library(mvtnorm)
library(VGAM)
library(invgamma)
library(Matrix)
Rcpp::sourceCpp("MCMC.cpp")
 
# load data
load("../Data/data.RData")
EE <- data$E
A <- matrix(0,data$p,data$p )
for (i in 1:nrow(EE)) {
  A[EE[i,1],EE[i,2]] <- 1
}

# set parameters 
a_sigma <- 1
b_sigma <- 1 
a_xi <- 2
b_xi <- 1
n_iter <- 1e4
burnin <- 1e3
k <- 200
uj1 <- .5
r <- 5
r_star <- 5
r_total <- r+r_star
u1_mean_prior <- u0_mean_prior <- matrix(0,r_total-2,1)
Sigma0 <- Sigma1 <- 0.001*diag(r_total-2)
a_q <- b_q <- 1
u00_mean_prior <- 0
sigma_u00 <- 0.1
u10_mean_prior <- 0
sigma_u10 <- 0.1
u00_star_mean_prior <- 0
sigma_u00_star <- 0.1
u10_star_mean_prior <- 0
sigma_u10_star <- 0.1
mu <- 2
nu <- 2

# run MCMC
fit <- MCMC(data$y,data$X,A,A_star,n_iter,mu,nu,a_sigma,b_sigma,a_xi,b_xi,
                  a_q,b_q,u00_mean_prior,sigma_u00,u10_mean_prior,sigma_u10,u00_star_mean_prior,sigma_u00_star,u10_star_mean_prior,sigma_u10_star,
                  r, r_star,uj1,u0_mean_prior,u1_mean_prior, Sigma0, Sigma1,k)
beta_postmed <- matrix(apply(fit$BETA[(burnin+1):n_iter,],2,median),ncol=1)
