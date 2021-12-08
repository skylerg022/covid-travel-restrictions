## Latent Variable Bayesian Probit Model for Ordinal Data
## Model 1: Homoskedasticity, no correlation
## Y_i = j if gamma_{j-1} <= z_i < gamma_j
## Z|beta   ~ N(x_i' \beta, 1)

library(tidyverse)
library(truncnorm)
library(cowplot)

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


# Gibb's sampler, uniform priors on cutoffs -------------------------------

n <- nrow(X)
J <- length(unique(control$int_travel_controls)) # Number of groups
ndraws <- 5000
draws <- matrix(0, nrow = ndraws, ncol = ncol(X) + J-2)
gamma <- c(-Inf, 0, 0.34, 0.8, 1.44, Inf)
gamma <- c(-Inf, 0, 1.1, 1.76, 2.7, Inf)
est_gamma <- 3:(length(gamma) - 1)
beta <- rep(0.5, ncol(X))

for (i in 1:ndraws) {
  z <- rtruncnorm(n, a = gamma[y], b = gamma[y+1], X%*%beta)
  sd_beta <- solve(t(X)%*%X)
  beta <- rmvnorm(1, sd_beta %*% t(X)%*%z, sd_beta) %>%
    t()
  for (j in est_gamma) {
    gamma[j] <- runif(1, max(max(z[z >= gamma[j-1] & z < gamma[j]]), gamma[j-1]),
                      min(min(z[z >= gamma[j] & z < gamma[j+1]]), gamma[j+1]))
  }
  draws[i,] <- c(beta, gamma[est_gamma])
}

coda::effectiveSize(draws)
plot(draws[,16], type = 'l')
plot(draws[,17], type = 'l')
plot(draws[,18], type = 'l')



# ROUND 2 -----------------------------------------------------------------

doit <- function(ndraws = 10000, warmup = 1000, eps = 0.1) {
  # Stuff for cutpoints
  # cuts <- c(-Inf, 0, # Q: Why the zero? Guess: Identifiability issues if estimating; only need to estimate three, not four
  #           cumsum(prop.table(table(y)))[-1] %>%
  #             qnorm(., mean = 1)) # I am just picking an initial value for cutpoints, right?
  # cuts <- c(-Inf, 0, 0.5, 1.1, 2.2, Inf)
  cuts <- c(-Inf, 0, 1, 1.5, 2.5, Inf) # Round 2 estimate
  # cuts <- c(-Inf, 0, 25, 50, 100, Inf)
  est.cuts <- 3:(length(cuts)-1)
  nc <- length(est.cuts)
  trans.cuts <- log( c(cuts[est.cuts[1]], diff(cuts[est.cuts])) )
  # eps <- c(0.2, 0.2, 0.5)
  # eps <- 0.1
  
  # prop.var <- matrix(c(eps, eps/3, eps/3,
  #                      eps/3, eps, eps/3,
  #                      eps/3, eps/3, eps),
  #                    nrow = 3, byrow = TRUE)
  b <- 0.04
  prop.var.chol <- matrix(c(sqrt(eps),  b,         b,
                            0,          sqrt(eps), b,
                            0,          0,         sqrt(eps)),
                          nrow = 3, byrow = TRUE)
  prop.var <- t(prop.var.chol) %*% prop.var.chol
  accept <- 0
  
  n <- nrow(X)
  J <- length(unique(control$int_travel_controls)) # Number of groups
  draws <- matrix(0, nrow = ndraws, ncol = ncol(X) + J-2)
  beta <- rnorm(ncol(X), 0, 2)
  
  for (i in 1:ndraws) {
    # Gibbs sample of z's and betas
    z <- rtruncnorm(n, a = cuts[y], b = cuts[y+1], X%*%beta)
    sd_beta <- solve(t(X)%*%X)
    beta <- rmvnorm(1, mean = sd_beta %*% t(X)%*%z, sigma = sd_beta) %>%
      t()
    
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
    
    
    # # SEPARATE UPDATES TO CUTPOINTS
    # # Metropolis algorithm for sampling cutpoints
    # prop.cuts <- cuts
    # prop.trans <- trans.cuts
    # for (j in 1:nc) {
    #   prop.trans[j] <- rnorm(1, sd=eps[j])
    #   prop.cuts[est.cuts] <- cumsum(sum(exp(prop.trans)))
    #   
    #   ## Evaluate the likelihood under proposed cutpoints
    #   mn <- X%*%beta
    #   # Proposed
    #   lwr <- prop.cuts[y]
    #   upr <- prop.cuts[y+1]
    #   prop.prob <- log(pnorm(upr, mn, 1) - pnorm(lwr, mn, 1)) # Q: What if I get -Inf?
    #   sum(prop.prob == -Inf)
    #   prop.prob[prop.prob == -Inf] <- -1e99
    #   # Current
    #   lwr <- cuts[y]
    #   upr <- cuts[y+1]
    #   cur.prob <- log(pnorm(upr, mn, 1) - pnorm(lwr, mn, 1))
    #   sum(cur.prob == -Inf)
    #   cur.prob[cur.prob == -Inf] <- -1e99
    #   llike.diff <- sum(prop.prob - cur.prob) #, na.rm = TRUE) # I am removing those observations that are infinity minus inf. Is this a problem?
    #   llike.diff
    #   ## Evaluate the Metropolis-Hasting Ratio
    #   MH <- llike.diff + sum(dnorm(prop.trans, 0, 10, log=TRUE) -
    #                            dnorm(trans.cuts, 0, 10, log=TRUE))
    #   ## Determine whether to accept or reject
    #   if (log(runif(1)) < MH) {
    #     cuts <- prop.cuts
    #     trans.cuts <- prop.trans
    #     accept[j] <- accept[j] + 1
    #   }
    # }
    
    draws[i,] <- c(beta, cuts[est.cuts])
  }
  
  cat('Acceptance rate:', accept/ndraws, '\n')
  
  draws <- draws[-(1:warmup),]
  colnames(draws) <- c(colnames(X), 
                       'delta3', 'delta4', 'delta5')
  return(draws)
}

# Run multiple chains
draws <- replicate(4, doit(ndraws = 20000, warmup = 2000, eps = 0.005),
                   simplify = FALSE) %>%
  do.call(rbind, .)

# Do thinning
idx <- seq(4, nrow(draws), by = 4)
draws_init <- draws
draws <- draws[idx,]

# Diagnostic Plots --------------------------------------------------------

coda::effectiveSize(draws)
acf(draws[,16])
acf(draws[,17])
acf(draws[,18])
plot(draws[,16], type = 'l')
plot(draws[,17], type = 'l')
plot(draws[,18], type = 'l')
apply(draws, 2, mean)


# Results -----------------------------------------------------------------

ci <- apply(draws, 2, function(x) quantile(x, c(0.025, 0.975))) %>%
  t() %>%
  round(3)
tbl_names <- colnames(draws) %>%
  str_replace('continent', 'Continent: ')
tbl_names[7:18] <- c('Reporting Full Vaccinations', 'Has New Cases',
                     'Log(Population Density)', 'Median Age',
                     'Log(GDP Per Capita)', 'Month (1 = Jan 2020)',
                     'Log(Average New Cases)', 'Proportion Fully Vaccinated',
                     'Time:Log(Average New Cases)', 'Delta3',
                     'Delta4', 'Delta5')
tbl_ci <- data.frame(coef = tbl_names,
                     est = colMeans(draws) %>%
                       round(3),
                     lower = ci[,1],
                     upper = ci[,2])
rownames(tbl_ci) <- NULL
write_csv(tbl_ci, 'data/tbl_ci.csv')

# Prediction: Taiwan in the future

control_tw <- control %>%
  filter(location == 'Taiwan')
p1 <- control_tw %>%
  ggplot(aes(time, avg_new_cases)) +
  geom_line() +
  labs(x = 'Time (1 = Jan 2020)',
       y = 'Average New Cases')
p2 <- control_tw %>%
  ggplot(aes(time, perc_fully_vaccinated)) +
  geom_line() +
  labs(x = 'Time (1 = Jan 2020)',
       y = 'Proportion Fully Vaccinated')
p3 <- control_tw %>%
  ggplot(aes(time, int_travel_controls)) +
  geom_line() +
  labs(x = 'Time (1 = Jan 2020)',
       y = 'Int. Travel Control Status') +
  scale_y_continuous(breaks= c(1, 2, 3, 4, 5))

# Plot
plot_grid(p1, p2, p3, nrow = 3)


# Prediction
tw_attrs <- control_tw[1,]
X_new <- c(1, 1, 0, 0, 0, 0, 1, 1, 
           tw_attrs$log_population_density,
           tw_attrs$median_age,
           tw_attrs$log_gdp_per_capita,
           28,
           2,
           0.8,
           28*2)
X_new[9:15] <- (X_new[9:15] - attr(X_scaled, "scaled:center")) / attr(X_scaled, "scaled:scale")


prob_open <- apply(draws, 1, function(pars) {
  beta <- pars[1:15]
  # The area under the curve less than this value is the probability
  #  that Taiwan's policies will be at most extreme screening
  delta3 <- pars[16] 
  muhat <- X_new %*% beta
  pnorm(delta3, muhat, 1)
})

mean(prob_open)
quantile(prob_open, c(0.025, 0.975))
ggplot() +
  geom_histogram(aes(prob_open, y = ..density..),
                 col = 'black',
                 fill = 'royalblue') +
  theme_bw() +
  labs(x = 'P(Y <= 1) (Screening or Less)',
       y = 'Density')





# Stan? -------------------------------------------------------------------

## Settings for Stan
nCores <- parallel::detectCores() - 1
options(mc.cores = nCores)          # Use all available cores.
rstan_options(auto_write = TRUE)    # Cache compiled code.
rstan_options(javascript = FALSE)     # FIXME: Effort to reduce RStudio Stan probs. May be helping.

# Fit Stan model
m <- stan_model(model_code = readLines("model5_simple.stan"))
data <- list(N = N,
             y = y, 
             X = X[,1, drop = FALSE],
             mu_0 = mu_0,
             tau2_0 = tau2_0)
fit <- sampling(m, data = data,
                iter = 500, warmup = 0,
                chains = 1) #nCores)
fit_extr <- extract(fit)

plot(fit_extr$beta[,1], type = 'l')
coda::effectiveSize(fit_extr$beta)

