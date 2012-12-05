\name{benchmark}
\Rdversion{1.1}
\alias{benchmark}
\title{
The Rsolnp Benchmark Problems Suite.
}
\description{
The function implements a set of benchmark problems against the MINOS solver of 
Murtagh and Saunders.
}
\usage{
benchmark(id = "Powell")
}
\arguments{
  \item{id}{
The name of the benchmark problem. A call to the function \code{\link{benchmarkids}}
will return the available benchmark problems.
}
}
\details{
The benchmarks were run on dual xeon server with 24GB of memory and windows 7 
operating system. The MINOS solver was used via the tomlab interface.
}
\value{
A data.frame containing the benchmark data. The description of the benchmark 
problem can be accessed throught the \code{description} attribute of 
the data.frame.
}
\references{
W.Hock and K.Schittkowski, \emph{Test Examples for Nonlinear Programming Codes}, 
Lecture Notes in Economics and Mathematical Systems. Springer Verlag, 1981.\cr 
Y.Ye, \emph{Interior algorithms for linear, quadratic, and linearly constrained 
non linear programming}, PhD Thesis, Department of EES Stanford University, 
Stanford CA.\cr
B.A.Murtagh and M.A.Saunders, \emph{MINOS 5.5 User's Guide, Report SOL 83-20R}, 
Systems Optimization Laboratory, Stanford University (revised July 1998).\cr
P. E. Gill, W. Murray, and M. A. Saunders, \emph{SNOPT An SQP algorithm for 
large-scale constrained optimization}, SIAM J. Optim., 12 (2002), pp.979-1006.
}
\author{
Alexios Ghalanos and Stefan Theussl\cr
Y.Ye (original matlab version of solnp)
}
\examples{
\dontrun{
benchmarkids()
benchmark(id = "Powell")
benchmark(id = "Alkylation")
benchmark(id = "Box")
benchmark(id = "RosenSuzuki")
benchmark(id = "Wright4")
benchmark(id = "Wright9")
benchmark(id = "Electron")
benchmark(id = "Permutation")
# accessing the description
test = benchmark(id = "Entropy")
attr(test, "description")
}
}
\keyword{optimize}