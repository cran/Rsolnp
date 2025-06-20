## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE--------------------------------------------------------
library(Rsolnp)

## -----------------------------------------------------------------------------
args(csolnp)

## -----------------------------------------------------------------------------
args(solnp)

## -----------------------------------------------------------------------------
prob <- solnp_problem_suite(number = 10)
sol <- csolnp(pars = prob$start, fn = prob$fn, gr = prob$gr, eq = prob$eq_fn, eq_b = prob$eq_b, 
              eq_jac = prob$eq_jac, ineq_fn = prob$ineq_fn, ineq_lower = prob$ineq_lower, 
              ineq_upper = prob$ineq_upper, ineq_jac = prob$ineq_jac, lower = prob$lower,
              upper = prob$upper)
print(prob$name)
print(c("convergence" = sol$convergence))
print(c("csolnp objective" = sol$objective, "best objective" = prob$best_fn))
print(sol$elapsed_time)

## -----------------------------------------------------------------------------
prob <- solnp_problem_suite(suite = "Other", number = 4)
sol <- csolnp(pars = prob$start, fn = prob$fn, gr = prob$gr, eq = prob$eq_fn, eq_b = prob$eq_b, 
              eq_jac = prob$eq_jac, ineq_fn = prob$ineq_fn, ineq_lower = prob$ineq_lower, 
              ineq_upper = prob$ineq_upper, ineq_jac = prob$ineq_jac, lower = prob$lower,
              upper = prob$upper)
print(prob$name)
print(c("convergence" = sol$convergence))
print(c("csolnp objective" = sol$objective, "best objective" = prob$best_fn))
print(sol$elapsed_time)

## ----warning=FALSE------------------------------------------------------------
prob <- solnp_problem_suite(number = 55)
sol <- csolnp_ms(prob$fn, gr = prob$gr, eq_fn = prob$eq_fn, eq_b = prob$eq_b, 
              eq_jac = prob$eq_jac, ineq_fn = prob$ineq_fn, ineq_lower = prob$ineq_lower, 
              ineq_upper = prob$ineq_upper, ineq_jac = prob$ineq_jac, lower = prob$lower,
              upper = prob$upper, n_candidates = 200, penalty = 1e5, 
              control = list(min_iter = 1000, max_iter = 100, tol = 1e-8),
              seed = 300)

print(paste0("Solution : ", round(sol$objective,3), " | Best Objective :", round(prob$best_fn,3)))
print(paste0("Equaility Violation : ", round(sol$kkt_diagnostics$eq_violation,3)))

