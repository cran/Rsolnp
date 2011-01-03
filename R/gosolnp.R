#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

#---------------------------------------------------------------------------------
# optimization by randomly restarted parameters using simulated parameter strategy
# Alexios Ghalanos 2010
#---------------------------------------------------------------------------------

# allowed distributions:
# 1: uniform (no confidence in the location of the parameter...somewhere in LB-UB space)
# 2: truncnorm (high confidence in the location of the parameter)
# 3: normal (Uncertainty in Lower and Upper bounds, but some idea about the dispersion about the location)
# ...

gosolnp = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, 
		ineqUB = NULL, LB = NULL, UB = NULL, control = list(), distr = rep(1, length(LB)), distr.opt = list(), 
		n.restarts = 1, n.sim = 20000, parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), rseed = NULL, ...)
{
	if( !is.null(pars) ) xnames = names(pars) else xnames = NULL
	if(is.null(control$trace)) trace = FALSE else trace = as.logical(control$trace)
	
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	
	# function requires both upper and lower bounds
	if(is.null(LB)) stop("\ngosolnp-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB)) stop("\ngosolnp-->error: the function requires upper parameter bounds\n", call. = FALSE)
	
	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars)) stop("\ngosolnp-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)
	
	np = 1:n
	
	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np)))) stop("\ngosolnp-->error: fixed indices out of bounds\n", call. = FALSE)
	}
	
	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean)) stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), 
						call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd)) stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), 
						call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean)) stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd)) stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	
	# initiate random search
	spars = .randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun, eqB = eqB,  
			ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, 
			distr = distr, distr.opt = distr.opt, n.restarts = n.restarts, n.sim = n.sim,
			trace = trace, rseed = rseed, xnames, parallel = parallel, parallel.control = parallel.control, ...)
	
	# initiate solver restarts
	if(trace) cat("\ngosolnp-->Starting Solver\n")
	solution = vector(mode = "list", length = n.restarts)
	if( parallel ){
		os = .Platform$OS.type
		if(is.null(parallel.control$pkg)){
			if( os == "windows" ) parallel.control$pkg = "snowfall" else parallel.control$pkg = "multicore"
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2
		} else{
			mtype = match(tolower(parallel.control$pkg[1]), c("multicore", "snowfall"))
			if(is.na(mtype)) stop("\nParallel Package type not recognized in parallel.control\n")
			parallel.control$pkg = tolower(parallel.control$pkg[1])
			if( os == "windows" && parallel.control$pkg == "multicore" ) stop("\nmulticore not supported on windows O/S\n")
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2 else parallel.control$cores = as.integer(parallel.control$cores[1])
		}
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			solution = mclapply(1:n.restarts, FUN = function(i) {
				xx =spars[i,]
				names(xx) = xnames
				ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun, eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, control = control, ...),
						silent = TRUE)
				if(inherits(ans, "try-error")){
					ans = list()
					ans$values = 1e10
					ans$convergence = 0
					ans$pars = rep(NA, length(xx))
				} 
				return( ans )
			}, mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			sfInit(parallel=TRUE, cpus = parallel.control$cores)
			sfExport("spars", "xnames", "fun", "eqfun", "eqB", "ineqfun","ineqLB", "ineqUB", "LB", "UB", "control", "...", local = TRUE)
			solution = sfLapply(as.list(1:n.restarts), fun = function(i) {
				xx = spars[i,]
				names(xx) = xnames
				ans = try(Rsolnp::solnp(pars = xx, fun = fun, eqfun = eqfun, eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, control = control, ...),
						silent = TRUE)
				if(inherits(ans, "try-error")){
					ans = list()
					ans$values = 1e10
					ans$convergence = 0
					ans$pars = rep(NA, length(xx))
				} 
				return( ans )
			})
			sfStop()
		}
	} else{
		solution = lapply(1:n.restarts, FUN = function(i){
		xx = spars[i,]
		names(xx) = xnames
		ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun, eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, control = control, ...),
				silent = TRUE)
		if(inherits(ans, "try-error")){
			ans = list()
			ans$values = 1e10
			ans$convergence = 0
			ans$pars = rep(NA, length(xx))
		} 
		return( ans )
		})
	}
	if(n.restarts>1){
		best = sapply(solution, FUN = function(x) if(x$convergence!=0) NA else x$values[length(x$values)])
		if(all(is.na(best))) stop("\ngosolnp-->Could not find a feasible starting point...exiting\n", call. = FALSE)
		nb = which(best == min(best, na.rm = TRUE))[1]
		solution = solution[[nb]]
		if(trace) cat("\ngosolnp-->Done!\n")
		solution$start.pars = spars[nb,]
		names(solution$start.pars) = xnames
		solution$rseed = rseed
	} else{
		solution = solution[[1]]
		solution$start.pars = spars[1,]
		names(solution$start.pars) = xnames
		solution$rseed = rseed
	}
	return(solution)
}

.randpars = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB, LB, UB, 
		distr, distr.opt, n.restarts, n.sim, trace = TRUE, rseed, xnames, parallel = FALSE,
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	if(trace) cat("\ngosolnp-->Calculating Random Initialization Parameters...")
	tmpenvir <- environment()
	assign("xnames", xnames, envir = tmpenvir)
	N = length(LB)
	
	rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	
	if(!is.null(fixed)) for(i in 1:length(fixed)) rndpars[,fixed[i]] = pars[fixed[i]]
	
	nf = 1:N
	
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	
	m = length(nf)
	
	set.seed(rseed)
	
	for(i in 1:m){
		j = nf[i]
		rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
				)
	}
	
	if(trace) cat("ok!\n")
	
	if(!is.null(ineqfun)){
		if(trace) cat("\ngosolnp-->Excluding Inequality Violations...\n")
		ineqv = matrix(NA, ncol = length(ineqLB), nrow = n.restarts*n.sim)
		
		# ineqv = t(apply(rndpars, 1, FUN = function(x) ineqfun(x)))
		if(length(ineqLB) == 1){
			ineqv = apply(rndpars, 1, FUN = function(x){
								names(x) = xnames
								ineqfun(x, ...)} )
			lbviol = sum(ineqv<ineqLB)
			ubviol = sum(ineqv>ineqUB)
			if( lbviol > 0 | ubviol > 0 ){
				vidx = c(which(ineqv<ineqLB), which(ineqv>ineqUB))
				vidx = unique(vidx)
				rndpars = rndpars[-c(vidx),]
				lvx = length(vidx)
			} else{
				vidx = 0
				lvx = 0
			}
		} else{
			ineqv = t(apply(rndpars, 1, FUN = function(x){
			names(x) = xnames
			ineqfun(x, ...)} ))
		
			# check lower and upper violations
			lbviol = apply(ineqv, 1, FUN = function(x) sum(any(x<ineqLB)))
			ubviol = apply(ineqv, 1, FUN = function(x) sum(any(x>ineqUB)))
			if( any(lbviol > 0) | any(ubviol > 0) ){
				vidx = c(which(lbviol>0), which(ubviol>0))
				vidx = unique(vidx)
				rndpars = rndpars[-c(vidx),]
				lvx = length(vidx)
				
			} else{
				vidx = 0
				lvx = 0
			}
		}
		if(trace) cat(paste("\n...Excluded ", lvx, "/",n.restarts*n.sim, " Random Sequences\n", sep = ""))
	}
	
	# evaluate function value
	if(trace) cat("\ngosolnp-->Evaluating Objective Function with Random Initialization Parameters...")
	if( parallel ){
		nx = dim(rndpars)[1]
		os = .Platform$OS.type
		if(is.null(parallel.control$pkg)){
			if( os == "windows" ) parallel.control$pkg = "snowfall" else parallel.control$pkg = "multicore"
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2
		} else{
			mtype = match(tolower(parallel.control$pkg[1]), c("multicore", "snowfall"))
			if(is.na(mtype)) stop("\nParallel Package type not recognized in parallel.control\n")
			parallel.control$pkg = tolower(parallel.control$pkg[1])
			if( os == "windows" && parallel.control$pkg == "multicore" ) stop("\nmulticore not supported on windows O/S\n")
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2 else parallel.control$cores = as.integer(parallel.control$cores[1])
		}
		if( parallel.control$pkg == "multicore" ){
			
			if(!exists("mclapply")){
				require('multicore')
			}
				evfun = mclapply(1:nx, FUN = function(i) .safefun(rndpars[i, ], fun, .env = tmpenvir, ...), 
						mc.cores = parallel.control$cores)
				evfun = as.numeric( unlist(evfun) )
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("rndpars", "fun", "tmpenvir", "...", local = TRUE)
			evfun = sfLapply(as.list(1:nx), fun = function(i) Rsolnp:::.safefun(rndpars[i, ], fun, .env = tmpenvir, ...))
			sfStop()
			evfun = as.numeric( unlist(evfun) )
		}
	} else{
		evfun = apply(rndpars, 1, FUN = function(x) .safefun(x, fun, .env = tmpenvir, ...))
	}
	if(trace) cat("ok!\n")
	
	if(trace) cat("\ngosolnp-->Sorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	ans = rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	if(trace) cat("ok!\n")
	
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if(trace){
		cat("\ngosolnp-->Starting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	rm(tmpenvir)
	return(ans)
}


.distr1 = function(LB, UB, n)
{
	runif(n, min = LB, max = UB)
}

.distr2 = function(LB, UB, n, mean, sd)
{
	rtruncnorm(n, a = as.double(LB), b = as.double(UB), mean = as.double(mean), sd = as.double(sd))
}

.distr3 = function(n, mean, sd)
{
	rnorm(n, mean = mean, sd = sd)
}