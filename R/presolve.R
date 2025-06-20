quick_stationary_feasible <- function(x0, lower, upper, ineq_fn, ineq_lower, ineq_upper, penalty = 1e4, maxit = 10) {
    standard_form_ineq_fn <- function(x) {
        h <- ineq_fn(x)
        c(ineq_lower - h, h - ineq_upper)
    }
    penalized <- function(x) {
        penalty_term <- sum(pmax(standard_form_ineq_fn(x), 0)^2)
        norm_term <- sum((x - x0)^2)
        norm_term + penalty * penalty_term
    }
    res <- optim(x0, penalized, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = maxit))
    res$par
}

sample_box_feasible <- function(lower, upper, eps = 1e-4) {
    effective_lower <- lower + eps * (upper > lower + 2*eps)
    effective_upper <- upper - eps * (upper > lower + 2*eps)
    is_too_narrow <- effective_lower >= effective_upper
    res <- numeric(length(lower))
    res[is_too_narrow] <- (lower[is_too_narrow] + upper[is_too_narrow]) / 2
    res[!is_too_narrow] <- effective_lower[!is_too_narrow] +
        (effective_upper[!is_too_narrow] - effective_lower[!is_too_narrow]) * runif(sum(!is_too_narrow))
    res
}


generate_feasible_starts <- function(n, fn = NULL, lower, upper, ineq_fn, ineq_lower, ineq_upper, maxit = 100, penalty = 1e4, eps = 1e-4, seed = NULL)
{
    if (!is.null(seed)) set.seed(seed)
    candidates <- vector("list", n)
    for (i in seq_len(n)) {
        x0 <- sample_box_feasible(lower, upper, eps = eps)
        if (!is.null(ineq_fn)) {
            feasible <- quick_stationary_feasible(x0, lower, upper, ineq_fn, ineq_lower, ineq_upper, maxit = maxit, penalty = penalty)
        } else {
            feasible <- x0
        }
        candidates[[i]] <- feasible
    }
    # Optionally, rank by objective value
    if (!is.null(fn)) {
        obj_vals <- sapply(candidates, fn)
        ord <- order(obj_vals, decreasing = FALSE)
        candidates <- candidates[ord]
    }
    out <- do.call(rbind, candidates)
    if (!is.null(names(lower))) colnames(out) <- names(lower)
    return(out)
}

