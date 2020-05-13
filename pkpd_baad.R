#!/usr/bin/env Rscript

library(rstan)
library(methods)
library(reshape2)
library(ggplot2)
library(plyr)
library(functional)
library(abind)
source("utils_baad.R")

options(rcpp.cache.dir="./rcpp-cache")
rstan_options(auto_write = TRUE)

theme_set(theme_bw() + theme(panel.grid.major=element_line(colour="grey80")))

set.seed(24234)

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

expose_stan_functions("pkpd_mvn_approxb.stan")

chains <- 4
iter <- 5000
warmup <- iter/2

# define true parameters
Lalpha_0 <- valogit(0.5)

Lalpha_s <- valogit(0.35)
lkappa   <- log(10/52)
lEmax    <- log(valogit(0.60) - Lalpha_s)

sigma_Lalpha_0 <- 0.2
sigma_lkappa   <- 0.5

sigma_y <- 5/100

## from the first group the first half is placebo, the rest is on
## treatment 1
J <- 100

J_prime <- 100

## let's say we cover a year and measure monthly
x <- seq(0, 52, length=13) / 52
T <- 13

## for simplicity assume the prime data set is at the same time-points
x_prime <- seq(0, 52, length=13) / 52
T_prime <- 13

## number of simulations per draw of the posterior to get
## approximated log-lik weight
J_tilde <- 500

## define weakly-informative prior
## phi is ordered as:
##  phi[1] = Lalpha_0;
##  phi[2] = Lalpha_s;
##  phi[3] = lkappa;
##  phi[4] = lEmax;
##  phi[5] = log(sigma_Lalpha_0);
##  phi[6] = log(sigma_lkappa);
##  phi[7] = log(sigma_y);
mu_phi_p <- c(0, -1, log(1/4), log(0.5), 0, 0, 0)
K_phi <- length(mu_phi_p)
Sigma_phi_p <- diag(c(2, 2, log(2), log(2), 1, 1, 1)^2)

K <- 1
delta <- array(0.1, dim=1)

mu_delta_p <- array(0, dim=1)
Sigma_delta_p <- matrix(1^2)

## simulated asymptotic values for each group
inv_valogit(Lalpha_s) ## placebo t=inf
inv_valogit(Lalpha_s + exp(lEmax)) ## drug=1 t=inf
inv_valogit(Lalpha_s + exp(lEmax + delta)) ## drug=2 t=inf


## simulate patient specific parameters
J_tot <- J + J_prime
subj_tot <- data.frame(id=1:J_tot
                      ,prime=rep(c(0,1), times=c(J, J_prime))
                      ,DRUG=0
                      ,eta_Lalpha_0=rnorm(J_tot, 0, 1)
                      ,eta_lkappa =rnorm(J_tot, 0, 1)
                   )

subj_tot$DRUG[(J/2+1):J] <- 1
subj_tot$DRUG[(J+1):J_tot] <- 2

subj       <- subset(subj_tot, prime==0)
subj_prime <- subset(subj_tot, prime==1)

pstream__ <- get_stream()
ysim       <- sim_posterior(subj, evaluate_model, TRUE)
ysim_prime <- sim_posterior(subj_prime, evaluate_model, TRUE)

dimnames(ysim) <- list(id=1:J, k=1:T)
dimnames(ysim_prime) <- list(id=(1+J):J_tot, k=1:T_prime)

M <- arrange(melt(ysim), id, k)
Mp <- arrange(melt(ysim_prime), id, k)
M$x <- x[M$k]
Mp$x<- x_prime[Mp$k]

M <- merge(rbind(M,Mp), subj_tot, by="id")

## look at the simulated means
plm <- ggplot(M, aes(x, value, group=id, colour=factor(DRUG))) +
    geom_line(alpha=0.3) +
        stat_summary(aes(group=DRUG), fun.data = "mean_cl_boot", position=position_dodge(width=0.02))


## create noise which we recycle for all of the different delta
## simulations
y_err       <- matrix(rnorm(J*T,             0, sigma_y), J,       T)
y_prime_err <- matrix(rnorm(J_prime*T_prime, 0, sigma_y), J_prime, T_prime)
y           <- ysim       + y_err

y_prime <- ysim_prime + y_prime_err

y_prime_bar <- colMeans(y_prime)

## plot the simulated data (not just the means)

M_obs <- arrange(melt(y), id, k)
Mp_obs <- arrange(melt(y_prime), id, k)
M_obs$x <- x[M_obs$k]
Mp_obs$x<- x_prime[Mp_obs$k]

M_obs <- merge(rbind(M_obs,Mp_obs), subj_tot, by="id")

color_scale <- scale_colour_brewer(type="qual", palette=2)

M_obs <- transform(M_obs,
                   Data=factor(prime, levels=c(0,1), labels=c("Internal", "External")),
                   Treatment=factor(DRUG, levels=0:2, labels=c("Placebo", "Treatment 1", "Treatment 2"))
                   )

M_obs_aggr <- ddply(M_obs, .(Data, Treatment, x), plyr::summarize, value=mean(value))

pl_obs <- ggplot(M_obs, aes(x, value, colour=Treatment)) +
    geom_line(alpha=0.4, aes(group=id)) + facet_grid(.~Data, labeller=label_both, drop=FALSE) +
        coord_cartesian(ylim=c(0.25,0.725)) + xlab("Time") +  ylab("Response") +
            ##stat_summary(aes(group=DRUG), fun.data = "mean_cl_boot", position=position_dodge(width=0.02)) +
                theme(legend.position=c(0.85, 0.25)) + color_scale

pl_obs

pl_reported <- pl_obs %+% subset(M_obs, prime==0) + geom_point(data=subset(M_obs_aggr, Data=="External"), size=I(4)) + geom_line(data=subset(M_obs_aggr, Data=="External"), size=I(1.5))

pl_reported

pl_local <- pl_obs %+% subset(M_obs, prime==0)

pl_local

ggsave("pkpd_simulated.png", pl_obs, width=7, height=3)
ggsave("pkpd_reported.png", pl_reported, width=7, height=3)
ggsave("pkpd_local.png", pl_local, width=7, height=3)

model <- stan_model("pkpd_mvn_approxb.stan")

## needed for new approach below
C <- 4 ## max # of chains
xi <- array(rnorm(C*2* 2*J_tilde), dim=c(C, 2, 2*J_tilde))

DRUG <- subj_tot$DRUG
DRUG_prime <- 2

## First profile the noisiness of the MVN approximation wrt to
## information about delta

S <- 51
delta_eval <- cbind(delta_1=seq(0,0.4,length=S))

R <- 1e3
delta_quants <- sim_delta(J_tilde=100, R, delta_eval)

delta_quants <- cbind(delta_quants, as.data.frame(delta_eval))

pl_delta_sim <- ggplot(delta_quants, aes(delta_1, q50)) +
    geom_ribbon(aes(ymin=q10, ymax=q90), fill="grey80") +
        geom_line(linetype=2, colour="blue") +
            ggtitle(paste("Approximate MVN Log-Likelihood")) +
                ylab(expression(paste("80% CI of log(p(", bar(y), "\'|", phi, "\'))"))) +
                    xlab(expression(delta[1]))
pl_delta_sim

ggsave("kpd_mvn_approx_delta_1b.pdf", pl_delta_sim, width=6, height=4)

## now simulate at delta[1] = 0 an increasing number of J_tilde to see
## scaling of 80% credible interval with J_tilde

Vsim_delta <- Vectorize(sim_delta, "J_tilde", SIMPLIFY=FALSE)

delta_disp <- matrix(c(0), 1)

J_tilde_eval <- round(10^seq(log10(50), log10(5e3), length=21))
delta_cov <- do.call(rbind, Vsim_delta(J_tilde_eval, R, delta_disp))

delta_cov <- transform(as.data.frame(delta_cov), J_tilde=J_tilde_eval, c80=q90-q10)

lfit <- lm(log10(c80) ~ 1 + offset(-0.5 * log10(J_tilde)), data=delta_cov)
a <- coef(lfit)[1]


pl_scaling <- ggplot(delta_cov, aes(J_tilde, c80)) + geom_line() +
    scale_x_log10(breaks=c(40, 100, 200, 400, 1000, 2000, 4000)) +
        scale_y_log10(breaks=c(1, 2, 4, 8, 16)) +
            geom_abline(intercept=a, slope=-0.5, linetype=2) +
            ggtitle(expression(paste("Scaling at ", delta[1], "=0 of the 80% CI"))) +
                ylab("80% CI Width") +
                    xlab(expression(tilde(J)))
pl_scaling

ggsave("kpd_mvn_approx_delta_1_scalingb.pdf", pl_scaling, width=6, height=4)


pl_inset <- ggplotGrob(pl_scaling + theme_bw(base_size=9))

ym <- -820
xm <- -0.015
pl_comb <- pl_delta_sim + annotation_custom(grob = pl_inset,
                                            xmin = xm, xmax = xm + 0.2,
                                            ymin = ym, ymax = ym + 500)
pl_comb


ggsave("mvn_approx_kpdb.pdf", pl_comb, width=6, height=4)


###
### Fit the model
###

## to the full data set first
base_data  <- list(J=J, T=T, K=K, K_phi=K_phi, y=y,
                   J_prime=J_prime, T_prime=T_prime, y_prime=y_prime,
                   x=x,
                   x_prime=x_prime,
                   mu_phi_p=mu_phi_p,
                   Sigma_phi_p=Sigma_phi_p,
                   mu_delta_p=mu_delta_p,
                   Sigma_delta_p=Sigma_delta_p,
                   J_tilde=J_tilde,
                   DRUG=DRUG,
                   C=C,
                   xi=xi)

fit_full <- sampling(model, data=modifyList(base_data, list(fit_all=1, fit_local=0)),
                     chains=4, iter=iter, warmup=warmup,
                     chain_id=1L)

parameters_ref <- c("Lalpha_0", "Lalpha_s", "lkappa", "lEmax", "sigma_Lalpha_0", "sigma_lkappa", "sigma_y", "delta")
parameters_local <- exclude(parameters_ref, "delta")

cat("FULL FIT:\n")
print(fit_full, pars=parameters_ref )

## Fit the model once just to local data
fit_local <- sampling(model, data=modifyList(base_data, list(fit_all=0, fit_local=1)),
                      chains=chains, iter=iter, warmup=warmup,
                      chain_id=1L)

cat("LOCAL FIT:\n")
print(fit_local, pars=parameters_ref)

## Fit the model to local data and the summary data
fit_mvn_approx <- sampling(model, data=modifyList(base_data, list(T=T, fit_all=0, fit_local=0)),
                           chains=chains, iter=iter, warmup=warmup,
                           chain_id=1L)

cat("MVN APPROX FIT:\n")
print(fit_mvn_approx, pars=parameters_ref)

sum_all <- list(full=summary(fit_full, pars=parameters_ref)$summary,
                local=summary(fit_local, pars=parameters_local)$summary,
                approx=summary(fit_mvn_approx, pars=parameters_ref)$summary)

sum_all <- ldply(sum_all, function(d) { vars <- rownames(d); transform(as.data.frame(d), var=vars)} )

sum_all <- transform(sum_all,
                     case=factor(.id, levels=c("local", "full", "approx")),
                     var=factor(as.character(var),
                         levels=rev(c("Lalpha_0", "Lalpha_s", "lkappa", "lEmax", "sigma_Lalpha_0", "sigma_lkappa", "sigma_y", "delta[1]") )))

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
    geom_pointrange(aes(ymin=X2.5., ymax=X97.5., colour=case, shape=case), position=position_dodge(width=0.8)) + coord_flip() + xlab(NULL) + ylab("Value [a.u.]") +
        geom_segment(aes(x=nvar-0.4, xend=nvar+0.4, y=value, yend=value), data=truth, alpha=I(0.8)) +
            scale_x_discrete(labels=rev(c(expression(L * alpha[0]), expression(L * alpha[s]), expression(l * kappa), expression(lEmax), expression(sigma[L * alpha[0] ]), expression(sigma[l * kappa ]), expression(sigma[y]), expression(delta[1]) ))) +
            scale_colour_brewer("Scenario", type="qual", palette="Dark2") +
        scale_shape_discrete("Scenario") +
            theme(legend.position=c(0.1, 0.175)) +
        ggtitle("Drug-disease model estimates, 95% CrI")
pl_est

ggsave("kpd_estb.pdf", pl_est, width=6, height=4)

sum_bias <- merge(truth, sum_all)

shade <- data.frame(l=0.5 + seq(0,7,by=2), m=1.5 + seq(0,7,by=2))

pl_bias <- ggplot(sum_bias) +
    geom_rect(data=shade, aes(xmin=l, xmax=m), fill="grey80", ymin=-Inf, ymax=Inf) +
            geom_hline(yintercept=0, alpha=I(0.7)) +
                geom_pointrange(aes(var, X50. - value, ymin=X2.5.- value, ymax=X97.5.- value, colour=case, shape=case), position=position_dodge(width=0.8)) + coord_flip() + xlab(NULL) + ylab("(Estimate - True Value) [a.u.]") +
                    scale_x_discrete(labels=rev(c(expression(L * alpha[0]), expression(L * alpha[s]), expression(l * kappa), expression(lEmax), expression(sigma[L * alpha[0] ]), expression(sigma[l * kappa ]), expression(sigma[y]), expression(delta[1]) ))) +
            scale_colour_brewer("Scenario", type="qual", palette="Dark2") +
        scale_shape_discrete("Scenario") +
           theme(legend.position=c(0.1, 0.2)) +
        ggtitle("Drug-disease model bias, 95% CrI")
pl_bias

ggsave("kpd_biasb.pdf", pl_bias, width=6, height=4)
ggsave("kpd_biasb.png", pl_bias, width=6, height=4)

ref <- subset(sum_all, case == "full")[,c("var", "X50.", "sd")]
names(ref) <- c("var", "value", "scale")

sum_rbias <- merge(ref, sum_all)

pl_rbias <- ggplot(subset(sum_rbias, case != "full"), aes(var, (X50. - value)/scale)) +
            geom_hline(yintercept=0, alpha=I(0.7)) +
                geom_pointrange(aes(ymin=(X2.5.- value)/scale, ymax=(X97.5.- value)/scale, colour=case, shape=case), position=position_dodge(width=0.8)) + coord_flip() + xlab(NULL) + ylab("Value [a.u.]") +
                    scale_x_discrete(labels=rev(c(expression(L * alpha[0]), expression(L * alpha[s]), expression(l * kappa), expression(lEmax), expression(sigma[L * alpha[0] ]), expression(sigma[l * kappa ]), expression(sigma[y]), expression(delta[1]) ))) +
            scale_colour_brewer("Scenario", type="qual", palette="Dark2") +
        scale_shape_discrete("Scenario") +
        ggtitle("Drug-disease model relative bias, 95% CrI")
pl_rbias

ggsave("kpd_rbiasb.pdf", pl_rbias, width=6, height=4)

sessionInfo()

save(list=ls(), file="pkpd-run.rda")
