## Latent Variable Bayesian Probit Model for Ordinal Data
## Model 1: Homoskedasticity, no correlation
## Y_i = j if gamma_{j-1} <= z_i < gamma_j
## Z|beta   ~ N(x_i' \beta, 1)

library(tidyverse)
library(truncnorm)

# set working directory if necessary
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

control <- read_csv('data/international_controls.csv') %>%
  mutate(log_population_density = log(population_density),
         log_gdp_per_capita = log(gdp_per_capita),
         has_new_cases = ifelse(avg_new_cases > 0, 1, 0),
         log_avg_new_cases = ifelse(has_new_cases, log(avg_new_cases), 0),
         has_new_deaths = ifelse(avg_new_deaths > 0, 1, 0),
         log_avg_new_deaths = ifelse(has_new_deaths, log(avg_new_deaths), 0))

X <- model.matrix(~ continent + vaccination_info + has_new_cases +
                    log_population_density + median_age + 
                    log_gdp_per_capita +
                    time*log_avg_new_cases + perc_fully_vaccinated,
                  data = control)
X_scaled <- scale(X[,9:15])
X[,9:15] <- X_scaled
y <- control$int_travel_controls + 1

# Gibb's sampler, uniform priors on cutoffs

n <- nrow(X)
J <- length(unique(control$int_travel_controls)) # Number of groups
ndraws <- 1000
draws <- matrix(0, nrow = ndraws, ncol = ncol(X) + J-2)
gamma <- c(-Inf, 0, 0.34, 0.8, 1.44, Inf)
est_gamma <- 3:(length(gamma) - 1)
beta <- rep(0.5, ncol(X))

for (i in 1:ndraws) {
  z <- rtruncnorm(n, a = gamma[y], b = gamma[y+1], X%*%beta)
  sd_beta <- solve(t(X)%*%X)
  beta <- rnorm(sd_beta %*% t(X)%*%z, sd_beta)
  for (j in est_gamma) {
    gamma[j] <- runif(1, max(max(z[z >= gamma[j-1] & z < gamma[j]]), gamma[j-1]),
                      min(min(z[z >= gamma[j] & z < gamma[j+1]]), gamma[j+1]))
  }
  draws[i,] <- c(beta, gamma[est_gamma])
}

coda::effectiveSize(draws)
plot(draws[,16], type = 'l')


# ROUND 2 -----------------------------------------------------------------

# Stuff for cutpoints
cuts <- c(-Inf, 0, # Q: Why the zero? Guess: Identifiability issues if estimating; only need to estimate three, not four
          cumsum(prop.table(table(y)))[-1] %>%
            qnorm(., mean = 1)) # I am just picking an initial value for cutpoints, right?
est.cuts <- 3:(length(cuts)-1)
nc <- length(est.cuts)
trans.cuts <- log( c(cuts[est.cuts[1]], diff(cuts[est.cuts])) )
# trans.cuts <- c(0, 2, 2.4)
eps <- 6
prop.var <- eps*diag(nc)
# prop.var <- diag(c(2, 4, 4))
# prop.var <- matrix(c(eps, -eps/3, -eps/3,
#                      -eps/3, eps, -eps/3,
#                      -eps/3, -eps/3, eps),
#                    nrow = 3, ncol = 3,
#                    byrow = TRUE)
accept <- 0

ndraws <- 1000
draws <- matrix(0, nrow = ndraws, ncol = ncol(X) + J-2)
# beta <- c( 1, rep(0.05, ncol(X)-1))
beta <- rnorm(ncol(X), 0, 2)

for (i in 1:ndraws) {
  # Gibbs sample of z's and betas
  z <- rtruncnorm(n, a = cuts[y], b = cuts[y+1], X%*%beta)
  sd_beta <- solve(t(X)%*%X)
  beta <- rnorm(sd_beta %*% t(X)%*%z, sd_beta)
  
  # Metropolis algorithm for sampling cutpoints
  prop.cuts <- cuts
  prop.trans <- trans.cuts + t(chol(prop.var))%*%rnorm(nc) #MVN Draw
  prop.cuts[est.cuts] <- cumsum(exp(prop.trans))
  
  ## Evaluate the likelihood under proposed cutpoints
  mn <- X%*%beta
  # Proposed
  lwr <- prop.cuts[y]
  upr <- prop.cuts[y+1]
  prop.prob <- log(pnorm(upr, mn, 1) - pnorm(lwr, mn, 1)) # Q: What if I get -Inf?
  sum(prop.prob == -Inf)
  prop.prob[prop.prob == -Inf] <- -1e99
  # Current
  lwr <- cuts[y]
  upr <- cuts[y+1]
  cur.prob <- log(pnorm(upr, mn, 1) - pnorm(lwr, mn, 1))
  sum(cur.prob == -Inf)
  cur.prob[cur.prob == -Inf] <- -1e99
  llike.diff <- sum(prop.prob - cur.prob) #, na.rm = TRUE) # I am removing those observations that are infinity minus inf. Is this a problem?
  llike.diff
  ## Evaluate the Metropolis-Hasting Ratio
  MH <- llike.diff + sum(dnorm(prop.trans, 0, 10, log=TRUE) -
                           dnorm(trans.cuts, 0, 10, log=TRUE))
  ## Determine whether to accept or reject
  if (log(runif(1)) < MH) {
    cuts <- prop.cuts
    trans.cuts <- prop.trans
    accept <- accept + 1
  }
  
  draws[i,] <- c(beta, cuts[est.cuts])
}

# Check diagnostic plots
coda::effectiveSize(draws)
accept
plot(draws[,16], type = 'l')
plot(draws[,17], type = 'l')
plot(draws[,18], type = 'l')
apply(draws, 2, mean)
