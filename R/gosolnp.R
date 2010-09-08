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
		n.restarts = 1, n.sim = 20000, use.multicore = FALSE, rseed = NULL, ...)
{
	if(is.null(tolower(control$trace))) trace = FALSE else trace = as.logical(control$trace)
	
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	
	# check mclapply
	if(use.multicore){
		if(!exists("mclapply")) stop("\ngosolnp-->error: to use multicore you must manually load the package\n", 
						call. = FALSE)
	}
	
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
			ineqfun = ineqfun, ineqLB - ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, 
			distr = distr, distr.opt = distr.opt, n.restarts = n.restarts, n.sim = n.sim,
			trace = trace, rseed = rseed,  ...)
	
	# initiate solver restarts
	if(trace) cat("\ngosolnp-->Starting Solver\n")
	solution = vector(mode = "list", length = n.restarts)
	if(use.multicore){
		solution = mclapply(1:n.restarts, FUN = function(i) solnp(pars = spars[i,], fun = fun, eqfun = eqfun, eqB = eqB, ineqfun = ineqfun,
							ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, control = control, ...))
	} else{
		solution = lapply(1:n.restarts, FUN = function(i) solnp(pars = spars[i,], fun = fun, eqfun = eqfun, eqB = eqB, ineqfun = ineqfun,
							ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB, control = control, ...))
	}
	if(n.restarts>1){
		best = sapply(solution, FUN = function(x) if(x$convergence!=0) NA else x$values[length(x$values)])
		if(all(is.na(best))) stop("\ngosolnp-->Could not find a feasible starting point...exiting\n", call. = FALSE)
		nb = which(best == min(best, na.rm = TRUE))[1]
		solution = solution[[nb]]
		if(trace) cat("\ngosolnp-->Done!\n")
		solution$start.pars = spars[nb,]
		solution$rseed = rseed
	} else{
		solution = solution[[1]]
		solution$start.pars = spars[1,]
		solution$rseed = rseed
	}
	return(solution)
}

.randpars = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB, LB, UB, 
		distr, distr.opt, n.restarts, n.sim, trace = TRUE, rseed, ...)
{
	if(trace) cat("\ngosolnp-->Calculating Random Initialization Parameters...")
	
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
		ineqv = t(apply(rndpars, 1, FUN = function(x) ineqfun(x, ...)))
		
		# check lower and upper violations
		lbviol = apply(ineqv, 1, FUN = function(x) sum(any(x<ineqLB)))
		ubviol = apply(ineqv, 1, FUN = function(x) sum(any(x>ineqUB)))
		vidx = c(which(lbviol>0), which(ubviol>0))
		vidx = unique(vidx)
		rndpars = rndpars[-c(vidx),]
		if(trace) cat(paste("\n...Excluded ", length(vidx), "/",n.restarts*n.sim, " Random Sequences\n", sep = ""))
	}
	
	# evaluate function value
	if(trace) cat("\ngosolnp-->Evaluating Objective Function with Random Initialization Parameters...")
	evfun = apply(rndpars, 1, FUN = function(x) .safefun(x, fun, ...))
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