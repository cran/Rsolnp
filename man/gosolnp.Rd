\name{gosolnp}
\alias{gosolnp}
\title{
Random Initialization and Multiple Restarts of the solnp solver.
}
\description{
When the objective function is non-smooth or has many local minima, it is hard to judge the
optimality of the solution, and this usually depends critically on the starting parameters.
This function enables the generation of a set of randomly chosen parameters from which to 
initialize multiple restarts of the solver (see note for details).
}
\usage{
gosolnp(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, 
LB = NULL, UB = NULL, control = list(), distr = rep(1, length(LB)), distr.opt = list(), n.restarts = 1, 
n.sim = 20000, parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), rseed = NULL, ...)
}
\arguments{
  \item{pars}{
The starting parameter vector. This is not required unless the fixed option is also used.
}
  \item{fixed}{
The numeric index which indicates those parameters which should stay fixed instead of being randomly generated.
}
  \item{fun}{
The main function which takes as first argument the parameter vector and returns
a single value.
}
  \item{eqfun}{
(Optional) The equality constraint function returning the vector of evaluated equality 
constraints.
}
  \item{eqB}{
(Optioal) The equality constraints.
}
  \item{ineqfun}{
(Optional) The inequality constraint function returning the vector of evaluated inequality 
constraints.
}
  \item{ineqLB}{
(Optional) The lower bound of the inequality constraints.
}
  \item{ineqUB}{
(Optional) The upper bound of the inequality constraints.
}
  \item{LB}{
The lower bound on the parameters. This is not optional in this function.
}
  \item{UB}{
The upper bound on the parameters. This is not optional in this function.
}
  \item{control}{
(Optional) The control list of optimization parameters. See \code{solnp} function documentation for details.
}
  \item{distr}{
A numeric vector of length equal to the number of parameters, indicating the choice of distribution 
to use for the random parameter  generation. Choices are uniform (1), truncated normal (2), and normal (3).
}
  \item{distr.opt}{
If any choice in \code{distr} was anything other than uniform (1), this is a list equal to the length of the
parameters with sub-components for the mean and sd, which are required in the truncated normal and
normal distributions.
}
  \item{n.restarts}{
The number of solver restarts required.
}
  \item{n.sim}{
The number of random parameters to generate for every restart of the solver. Note that there will always
be significant rejections if inequality bounds are present. Also, this choice should also be motivated
by the width of the upper and lower bounds.
}
  \item{parallel}{ 
Whether to make use of parallel processing on multicore systems.
}
  \item{parallel.control}{ 
The parallel control options including the type of package for performing the parallel calculations 
(\sQuote{multicore} for non-windows O/S and \sQuote{snowfall} for all O/S), and the number of cores to make use of.
}
  \item{rseed}{
(Optional) A seed to initiate the random number generator, else system time will be used.
}
  \item{\dots}{
(Optional) Additional parameters passed to the main, equality or inequality functions
}
}
\details{
Given a set of lower and upper bounds, the function generates, for those parameters not set as fixed,
random values from one of the 3 chosen distributions. The parameters are then passed to the inequality
function if that is present, and any parameters violating the upper or lower inequality bounds are
excluded. The parameters are then passed to the objective function, which is evaluated at those
parameters, the resulting values sorted, and the best N (N = random.restart) parameter vectors 
(corresponding to the best N objective function values) chosen in order to initialize the solver.
}
\value{
A list containing the following values:
\item{pars}{Optimal Parameters.}
\item{convergence }{Indicates whether the solver has converged (0) or not (1).}
\item{values}{Vector of function values during optimization with last one the
value at the optimal.}
\item{lagrange}{The vector of Lagrange multipliers.}
\item{hessian}{The Hessian at the optimal solution.}
\item{ineqx0}{The estimated optimal inequality vector of slack variables used for 
transforming the inequality into an equality constraint.}
\item{nfuneval}{The number of function evaluations.}
\item{elapsed}{Time taken to compute solution.}
\item{start.pars}{The parameter vector used to start the solver}
}
\references{
Y.Ye, \emph{Interior algorithms for linear, quadratic, and linearly constrained non 
linear programming}, PhD Thesis, Department of EES Stanford University, Stanford CA.\cr
Hu, X. and Shonkwiler, R. and Spruill, M.C. \emph{Random Restarts in Global Optimization},
1994, Georgia Institute of technology, Atlanta.
}
\author{
Alexios Ghalanos and Stefan Theussl\cr
Y.Ye (original matlab version of solnp)
}
\note{
The choice of which distribution to use for randomly sampling the parameter space should be driven by the
user's knowledge of the problem and confidence or lack thereof of the parameter distribution. The uniform distribution
indicates a lack of confidence in the location or dispersion of the parameter, while the truncated normal indicates
a more confident choice in both the location and dispersion. On the other hand, the normal indicates perhaps a lack
of knowledge in the upper or lower bounds, but some confidence in the location and dispersion of the parameter.
In using choices (2) and (3) for \code{distr}, the \code{distr.opt} list must be supplied with \code{mean} and 
\code{sd} as subcomponents for those parameters not using the uniform (the examples section hopefully clarifies the usage).
}
\examples{
\dontrun{
# Distributions of Electrons on a Sphere Problem
# Given n electrons, find the equilibrium state distribution (of minimal Coulomb potential) 
# of the electrons positioned on a conducting sphere. This model is from the COPS benchmarking suite.
# See http://www-unix.mcs.anl.gov/~more/cops/.	
	gofn = function(dat, n)
	{
		
		x = dat[1:n]
		y = dat[(n+1):(2*n)]
		z = dat[(2*n+1):(3*n)]
		ii = matrix(1:n, ncol = n, nrow = n, byrow = TRUE)
		jj = matrix(1:n, ncol = n, nrow = n)
		ij = which(ii<jj, arr.ind = TRUE)
		i = ij[,1]
		j = ij[,2]
		#  Coulomb potential
		potential = sum(1.0/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2))
		potential
	}
	
	goeqfn = function(dat, n)
	{
		x = dat[1:n]
		y = dat[(n+1):(2*n)]
		z = dat[(2*n+1):(3*n)]
		apply(cbind(x^2, y^2, z^2), 1, "sum")
	}
	
	n = 25
	LB = rep(-1, 3*n)
	UB = rep(1,  3*n)
	eqB = rep(1, n)
	ans = gosolnp(pars  = NULL, fixed = NULL, fun = gofn, eqfun = goeqfn, eqB = eqB, LB = LB, UB = UB, 
			control = list(outer.iter = 100, trace = 1), distr = rep(1, length(LB)), distr.opt = list(), 
			n.restarts = 2, n.sim = 20000, rseed = 443, n = 25)
	# should get a function value around 243.813
}
}
\keyword{optimize}