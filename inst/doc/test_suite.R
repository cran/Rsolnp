## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE--------------------------------------------------------
library(Rsolnp)

## -----------------------------------------------------------------------------
head(solnp_problems_table())

## -----------------------------------------------------------------------------
prob <- solnp_problem_suite(suite = "Hock-Schittkowski", number = 60)
sol <- csolnp(prob$start, fn = prob$fn, gr = prob$gr, eq_fn = prob$eq_fn, 
              eq_b = prob$eq_b, eq_jac = prob$eq_jac, ineq_fn = prob$ineq_fn, 
              ineq_lower = prob$ineq_lower, ineq_upper = prob$ineq_upper, 
              ineq_jac = prob$ineq_jac, lower = prob$lower, upper = prob$upper, 
              control = list(trace = 0, min_iter = 1000, max_iter = 300, tol = 1e-8, rho = 1))
print(c("solnp" = sol$objective, "benchmark" = prob$best_fn))
print(rbind("solnp" = sol$pars, "benchmark" = prob$best_par))

