\name{solnp}
\alias{solnp}
\title{
Nonlinear optimization using augmented lagrange method.
}
\description{
The solnp function is based on the solver by Yinyu Ye which solves the general nonlinear
programming problem:\cr
\deqn{\min  f(x)}{min  f(x)}
\deqn{\mathrm{s.t.}}{s.t.}
\deqn{g(x) = 0}{g(x) = 0}
\deqn{l_h \leq h(x) \leq u_h}{l[h] <= h(x) <= u[h]}
\deqn{l_x \leq x \leq u_x}{l[x] <= x <= u[x]}
where, \eqn{f(x)}, \eqn{g(x)} and \eqn{h(x)} are smooth functions.
}
\usage{
solnp(pars, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, 
LB = NULL, UB = NULL, control = list(), ...)
}
\arguments{
  \item{pars}{
The starting parameter vector.
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
(Optional) The equality constraints.
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
(Optional) The lower bound on the parameters.
}
  \item{UB}{
(Optional) The upper bound on the parameters.
}
  \item{control}{
(Optional) The control list of optimization parameters. See below for details.
}
  \item{\dots}{
(Optional) Additional parameters passed to the main, equality or inequality functions. Note
that the main and constraint functions must take the exact same arguments, irrespective of
whether they are used by all of them.
}
}
\details{
The solver belongs to the class of indirect solvers and implements the augmented 
Lagrange multiplier method with an SQP interior algorithm.
}
\value{
A list containing the following values:
\item{pars}{Optimal Parameters.}
\item{convergence }{Indicates whether the solver has converged (0) or not (1 or 2).}
\item{values}{Vector of function values during optimization with last one the
value at the optimal.}
\item{lagrange}{The vector of Lagrange multipliers.}
\item{hessian}{The Hessian of the augmented problem at the optimal solution.}
\item{ineqx0}{The estimated optimal inequality vector of slack variables used for 
transforming the inequality into an equality constraint.}
\item{nfuneval}{The number of function evaluations.}
\item{elapsed}{Time taken to compute solution.}
}
\section{Control}{
\describe{
\item{rho}{This is used as a penalty weighting scaler for infeasibility in the augmented
objective function. The higher its value the more the weighting to bring the solution into the 
feasible region (default 1). However, very high values might lead to numerical ill conditioning 
or significantly slow down convergence.}
\item{outer.iter}{Maximum number of major (outer) iterations (default 400).}
\item{inner.iter}{Maximum number of minor (inner) iterations (default 800).}
\item{delta}{Relative step size in forward difference evaluation (default 1.0e-7).}
\item{tol}{ Relative tolerance on feasibility and optimality (default 1e-8).}
\item{trace}{The value of the objective function and the parameters is printed at every
major iteration (default 1).}
}}
\references{
Y.Ye, \emph{Interior algorithms for linear, quadratic, and linearly constrained non 
linear programming}, PhD Thesis, Department of EES Stanford University, Stanford CA.
}
\author{
Alexios Ghalanos and Stefan Theussl\cr
Y.Ye (original matlab version of solnp)
}
\note{
The control parameters \code{tol} and \code{delta} are key in getting any possibility of 
successful convergence, therefore it is suggested that the user change these appropriately 
to reflect their problem specification.\cr
The solver is a local solver, therefore for problems with rough surfaces and many local
minima there is absolutely no reason to expect anything other than a local solution.
}
\examples{
# From the original paper by Y.Ye
# see the unit tests for more....
#---------------------------------------------------------------------------------
# POWELL Problem
fn1=function(x)
{
	exp(x[1]*x[2]*x[3]*x[4]*x[5])
}

eqn1=function(x){
	z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
	z2=x[2]*x[3]-5*x[4]*x[5]
	z3=x[1]*x[1]*x[1]+x[2]*x[2]*x[2]
	return(c(z1,z2,z3))
}


x0 = c(-2, 2, 2, -1, -1)
powell=solnp(x0, fun = fn1, eqfun = eqn1, eqB = c(10, 0, -1))
}
\keyword{optimize}