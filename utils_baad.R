## numerically stable summation of logs on the natural scale
log_sum_exp <- function(x){ 
   xmax <- which.max(x) 
   log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax] 
} 

## exclude elements from a list
exclude <- function(set, elems, keys=set) {
    set[-match(elems, keys)]
}

## extract from a posterior given as as list a specifc draw. Assumes
## that the first dimension of each list entry is the iteration.
extract_draw <- function(sims, draw) lapply(sims, asub, idx=draw, dim=1)

## applies a function over each entry of the posterior if
## vectorized=FALSE; for vectorized=TRUE the function is assumed to
## perform the simulation in a single sweep. Note that all arguments
## to the function are automatically deduced from it's formals and
## that all arguments which are not in the sims list are searched in
## the global environment.
sim_posterior <- function(sims, fun, vectorized=FALSE, res_type, envir) {
    args <- setdiff(names(formals(fun)), "seed")

    from_draw <- intersect(args, names(sims))
    from_env  <- setdiff(args, names(sims))

    if (missing(envir))
        envir <- parent.frame()

    sims <- sims[from_draw]
    aux <- mget(from_env, envir=envir)

    if(!vectorized) {
        S <- NROW(sims[[1]])
        calc_draw <- function(i) do.call(fun, c(aux, extract_draw(sims, i)))
        if(missing(res_type))
            res_type <- calc_draw(1)
        return(t(vapply(1:S, calc_draw, res_type)))
    } else {
        return(do.call(fun, c(sims, aux)))
    }
}


relabel <- function(x) match(x, unique(x))


## simulates for different delta's the mvn approximation which needs
## to be defined as mvn_approx_log in the global namespace
sim_delta <- function(J_tilde, R, delta, probs=c(0.1, 0.5, 0.9)) {
    ##J_tilde number of fake patients per combination
    S <- nrow(delta)  ## different number of deltas
    delta_eval <- delta[rep(1:S, each=R),,drop=FALSE]
    ## R is the number of replications per delta
    sim_delta <- list(
        delta=delta_eval
       ,xi=array(rnorm(R*S*2*J_tilde), dim=c(R*S,2,J_tilde))
        )

    mvn_delta_eval <- sim_posterior(sim_delta, mvn_approx_lpdf, envir=globalenv())

    delta_sim <- matrix(mvn_delta_eval, nrow=R, ncol=S)

    delta_quants <- t(apply(delta_sim, 2, quantile, probs=probs))

    colnames(delta_quants) <- paste0("q", round(probs * 100))
    as.data.frame(delta_quants)
}

amean <- function(x) {
    nd <- length(dim(x))
    if(nd == 1)
        return(mean(x))
    if(nd == 2)
        return(colMeans(x))
    stop("Unkown object")
}
