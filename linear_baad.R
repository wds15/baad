#!/usr/bin/env Rscript

library(rstan)
library(methods)
library(abind)
library(mvtnorm)
library(plyr)
library(ggplot2)
source("utils_baad.R")

options(rcpp.cache.dir="./rcpp-cache")
rstan_options(auto_write = TRUE)

theme_set(theme_bw() + theme(panel.grid.major=element_line(colour="grey80")))

set.seed(242345)

## choose whatever is appropiate for your computing environment
## useful on laptop
##cores <- detectCores()
## needed on cluster
cores <- as.numeric(Sys.getenv("NSLOTS"))

if(is.na(cores))
    cores <- 4

options(mc.cores = cores)

## programs must be run with more than one chain to make the automatic
## definition of CHAIN_ID work when calling stan
if(cores == 1)
    stop("Program must run with more than 1 chain!")

cat("Using", cores, "CPU cores\n")

expose_stan_functions("linear_mvn_approx.stan")

chains <- 4
iter <- 5000
warmup <- iter/2

# Define the hyperparameters of the model
mu_a <- c(.5, -.2)
sigma_a <- c(.1, .1)
Sigma_a <- diag(sigma_a^2)
sigma_y <- .05
K <- length(mu_a)
beta <- -.1

# Set up the experimental conditions
J <- 100
T <- 13
x <- (0:(T-1)) /(T-1)

# Simulate the data
a <- rmvnorm(J, mu_a, Sigma_a)
y <- array(NA, c(J,T))
for (k in 1:T)
  y[,k] <- rnorm(J, a[,1] + a[,2]*x[k] + beta*x[k]^2, sigma_y)

# Set up the conditions for the external data
J_prime <- 100
T_prime <- 13
x_prime <- (0:(T_prime-1)) / (T_prime-1)

# Simulate external data and compute averages
delta <- c(.1, .1)
a_prime <- rmvnorm(J_prime, mu_a + delta, Sigma_a)
y_prime <- array(NA, c(J_prime,T_prime))
for (k in 1:T_prime)
  y_prime[,k] <- rnorm(J_prime, a_prime[,1] + a_prime[,2]*x_prime[k] + beta*x_prime[k]^2, sigma_y)
y_prime_bar <- colMeans(y_prime)
true <- c(mu_a, beta, sigma_a, sigma_y, delta)

# Set noninformative priors
K_phi <- 6
mu_phi_p <- rep(0,K_phi)
Sigma_phi_p <- diag(K_phi)
mu_delta_p <- rep(0,K)
Sigma_delta_p <- diag(K)

# number of fake patients for simulations
J_tilde <- 500

parameters_ref <- c("mu_a", "beta", "sigma_a", "sigma_y", "delta")

## FIRST EXPLORE HOW MANY FAKE SIMULATIONS WE NEED

dens_correct <- function(y_prime_bar, mu_a, delta, b) {
    a <- mu_a + delta
    X_prime <- cbind(1, x_prime)
    y_correct_var <- (X_prime %*% diag(sigma_a^2) %*% t(X_prime)  + diag(sigma_y^2, T_prime))/J_prime
    y_correct_bar <- a[1] + a[2] * x_prime + b*x_prime^2
    dmvnorm(y_prime_bar, mean=y_correct_bar, sigma=y_correct_var, log=TRUE)
}

dens_correct_arg <- function(y_prime_bar, mu_a, delta_1, delta_2, b) {
    return(dens_correct(y_prime_bar, mu_a, c(delta_1, delta_2), b))
}

Vdens_correct_arg <- Vectorize(dens_correct_arg, c("delta_1", "delta_2"))

S <- 51
delta_eval <- cbind(delta_1=delta[1], delta_2=seq(0,0.22,length=S))

R <- 1e3
pstream__ <- get_stream()
delta_quants <- sim_delta(J_tilde=100, R, delta_eval)
delta_quants <- cbind(delta_quants, as.data.frame(delta_eval))

pl_delta_sim <- ggplot(delta_quants, aes(delta_2, q50)) +
    geom_ribbon(aes(ymin=q10, ymax=q90), fill="grey80") +
        geom_line(linetype=2, colour="blue") +
        stat_function(fun=function(x) Vdens_correct_arg(y_prime_bar, mu_a, delta[1], x, beta) ) +
            ggtitle(paste("Approximate MVN Log-Likelihood")) +
                ylab(expression(paste("80% CI of log(p(", bar(y), "\'|", phi, "\'))"))) +
                    xlab(expression(delta[2]))
pl_delta_sim

ggsave("mvn_approx_delta_2.pdf", pl_delta_sim, width=6, height=4)

## now simulate at delta[2] = 0 an increasing number of J_tilde to see
## scaling of 80% confidence interval with J_tilde

Vsim_delta <- Vectorize(sim_delta, "J_tilde", SIMPLIFY=FALSE)

delta_disp <- matrix(c(delta[1], 0), 1)

J_tilde_eval <- round(10^seq(log10(50), log10(5e3), length=21))
delta_cov <- do.call(rbind, Vsim_delta(J_tilde_eval, R, delta_disp))

delta_cov <- transform(as.data.frame(delta_cov), J_tilde=J_tilde_eval, c80=q90-q10)

lfit <- lm(log10(c80) ~ 1 + offset(-0.5 * log10(J_tilde)), data=delta_cov)
a <- coef(lfit)[1]

pl_scaling <- ggplot(delta_cov, aes(J_tilde, c80)) + geom_line() +
    scale_x_log10(breaks=c(40, 100, 200, 400, 1000, 2000, 4000)) +
        scale_y_log10(breaks=c(1, 2, 4, 10, 20, 40)) +
            geom_abline(intercept=a, slope=-0.5, linetype=2) +
            ggtitle(expression(paste("Scaling at ", delta[2], "=0 of the 80% CI"))) +
                ylab("80% CI Width") +
                    xlab(expression(tilde(J)))
pl_scaling

ggsave("mvn_approx_delta_2_scaling.pdf", pl_scaling, width=6, height=4)


pl_inset <- ggplotGrob(pl_scaling + theme_bw(base_size=9))

ym <- -50
xm <- -0.01
pl_comb <- pl_delta_sim + annotation_custom(grob = pl_inset,
                                            xmin = xm, xmax = xm + 0.11,
                                            ymin = ym, ymax = ym + 52.5)
pl_comb

ggsave("mvn_approx_linear.pdf", pl_comb, width=6, height=4)


## NOW FIT THE MODELS
model <- stan_model("linear_mvn_approx.stan")

## needed for new approach below
C <- 4 ## max # of chains
xi <- array(rnorm(C*K* 2*J_tilde), dim=c(C, K, 2*J_tilde))

## Fit the model to the complete data
base_data  <- list(J=J, T=T, K=K, K_phi=K_phi, y=y,
                   J_prime=J_prime, T_prime=T_prime, y_prime=y_prime,
                   x=x,
                   x_prime=x_prime,
                   mu_phi_p=mu_phi_p,
                   Sigma_phi_p=Sigma_phi_p,
                   mu_delta_p=mu_delta_p,
                   Sigma_delta_p=Sigma_delta_p,
                   J_tilde=J_tilde,
                   C=C,
                   xi=xi)

fit_full <- rstan::sampling(model, data=modifyList(base_data, list(fit_all=1, fit_local=0)),
                     chains=chains, iter=iter, warmup=warmup,
                     chain_id=1L)

cat("FULL FIT:\n")
print(fit_full, pars=parameters_ref )

## Fit the model to y and y_prime_bar
fit_correct <- sampling(model, data=modifyList(base_data, list(fit_all=-1, fit_local=0)),
                        chains=chains, iter=iter, warmup=warmup,
                        chain_id=1L)

cat("INTEGRATED FIT:\n")
print(fit_correct, pars=parameters_ref)

## Fit the model once just to local data
fit_local <- sampling(model, data=modifyList(base_data, list(fit_all=0, fit_local=1)),
                      chains=chains, iter=iter, warmup=warmup,
                      chain_id=1L)

cat("LOCAL FIT:\n")
print(fit_local, pars=parameters_ref)

## Fit the model to local data and the summary data
fit_mvn_approx <- sampling(model, data=modifyList(base_data, list(fit_all=0, fit_local=0)),
                           chains=chains, iter=iter, warmup=warmup,
                           chain_id=1L)

cat("MVN APPROX FIT:\n")
print(fit_mvn_approx, pars=parameters_ref)



sum_all <- list(full=summary(fit_full, pars=parameters_ref)$summary,
                local=summary(fit_local, pars=exclude(parameters_ref, "delta"))$summary,
                integrated=summary(fit_correct, pars=parameters_ref)$summary,
                approx=summary(fit_mvn_approx, pars=parameters_ref)$summary)


sum_all <- ldply(sum_all, function(d) { vars <- rownames(d); transform(as.data.frame(d), var=vars)} )

sum_all <- transform(sum_all,
                     case= factor(.id, levels=c("local", "full", "approx", "integrated")),
                     var=factor(as.character(var),
                         levels=rev(c("mu_a[1]", "mu_a[2]", "beta", "sigma_a[1]", "sigma_a[2]", "sigma_y", "delta[1]", "delta[2]") )))


## calculate summary of Neff + sem to check MCMC precision is great enough
ddply(sum_all, .(.id), summarize,
      neff_gmean=exp(mean(log(n_eff))),
      neff_min=min(n_eff),
      neff_minvar=var[which.min(n_eff)],
      sem_gmean=exp(mean(log(se_mean))),
      sem_max=max(se_mean),
      sem_maxvar=var[which.max(se_mean)]
      )

sum_all[,c(".id", "var", "n_eff", "se_mean")]

nf <- c(1,5:9)
## extract maximal sem:
cat("Maximal SEM:\n")
m <- which.max(sum_all$se_mean)
sum_all[m,-nf]

cat("Maximal SEM without delta:\n")
sum_all_sub <- sum_all[grep("delta", sum_all$var, invert=TRUE),]
m <- which.max(sum_all_sub$se_mean)
sum_all_sub[m,-nf]

truth <- data.frame(var=subset(sum_all, case == "full")$var, value=NA)

truth$value <- sapply(as.character(truth$var), function(v) eval(parse(text=v)))
truth$nvar <- as.numeric(truth$var)

pl_est <- ggplot(sum_all, aes(var, X50.)) +
    geom_pointrange(aes(ymin=X2.5., ymax=X97.5., colour=case, shape=case), position=position_dodge(width=0.8)) + coord_flip(ylim=c(-0.25,0.55)) + xlab(NULL) + ylab("Value [a.u.]") +
        scale_x_discrete(labels=rev(c(expression(mu[alpha*1]), expression(mu[alpha*2]), expression(beta), expression(sigma[alpha*1]), expression(sigma[alpha*2]), expression(sigma[y]), expression(delta[1]), expression(delta[2]))), drop=FALSE) +
        geom_segment(aes(x=nvar-0.4, xend=nvar+0.4, y=value, yend=value), data=truth, alpha=I(0.8)) +
            scale_colour_brewer("Scenario", type="qual", palette="Dark2") +
        scale_shape_discrete("Scenario") +
            theme(legend.position=c(0.1, 0.2)) +
        ggtitle("Linear model estimates, 95% CrI")
pl_est

ggsave("linear_est.pdf", pl_est, width=6, height=4)

sum_bias <- merge(truth, sum_all)

shade <- data.frame(l=0.5 + seq(0,7,by=2), m=1.5 + seq(0,7,by=2))

pl_bias <- ggplot(sum_bias) +
    geom_rect(data=shade, aes(xmin=l, xmax=m), fill="grey80", ymin=-Inf, ymax=Inf) +
            geom_hline(yintercept=0, alpha=I(0.7)) +
    geom_pointrange(aes(var, X50. - value, ymin=X2.5.- value, ymax=X97.5.- value, colour=case, shape=case), position=position_dodge(width=0.8)) + coord_flip(ylim=c(-0.075, 0.075)) + xlab(NULL) + ylab("(Estimate - True Value) [a.u.]") +
        scale_x_discrete(labels=rev(c(expression(mu[alpha*1]), expression(mu[alpha*2]), expression(beta), expression(sigma[alpha*1]), expression(sigma[alpha*2]), expression(sigma[y]), expression(delta[1]), expression(delta[2]))), drop=FALSE) +
            scale_colour_brewer("Scenario", type="qual", palette="Dark2") +
        scale_shape_discrete("Scenario") +
            theme(legend.position=c(0.1, 0.21)) +
        ggtitle("Linear model bias, 95% CrI")
pl_bias

ggsave("linear_bias.pdf", pl_bias, width=6, height=4)

ref <- subset(sum_all, case == "full")[,c("var", "X50.", "sd")]
names(ref) <- c("var", "value", "scale")

sum_rbias <- merge(ref, sum_all)

pl_rbias <- ggplot(subset(sum_rbias, case != "full"), aes(var, (X50. - value)/scale)) +
            geom_hline(yintercept=0, alpha=I(0.7)) +
    geom_pointrange(aes(ymin=(X2.5.- value)/scale, ymax=(X97.5.- value)/scale, colour=case, shape=case), position=position_dodge(width=0.8)) + coord_flip() + xlab(NULL) + ylab("(Estimate - Estimate_full)/scale_full [a.u.]") +
        scale_x_discrete(labels=rev(c(expression(mu[alpha*1]), expression(mu[alpha*2]), expression(beta), expression(sigma[alpha*1]), expression(sigma[alpha*2]), expression(sigma[y]), expression(delta[1]), expression(delta[2]))), drop=FALSE) +
            scale_colour_brewer("Scenario", type="qual", palette="Dark2") +
        scale_shape_discrete("Scenario") +
        ggtitle("Linear model estimates - relative bias, 95% CrI")
pl_rbias

ggsave("linear_rbias.pdf", pl_rbias, width=6, height=4)

sessionInfo()

save(list=ls(), file="linear-run.rda")
