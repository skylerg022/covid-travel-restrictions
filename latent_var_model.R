## Latent Variable Bayesian Probit Model for Ordinal Data
## Model 1: Homoskedasticity, no correlation
## Y_i = j if gamma_{j-1} <= z_i < gamma_j
## Z|beta   ~ N(x_i' \beta, 1)

library(tidyverse)
library(truncnorm)
library(cowplot)
library(mvtnorm)

set.seed(187)

# set working directory if necessary
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Set ggplot theme
old_theme <- theme_set(theme_bw())
# Picture saving settings
pic_width <- 7
pic_height <- 4
pic_unit <- 'in'


# Read in data
control <- read_csv('data/international_controls.csv') %>%
  mutate(log_population_density = log(population_density),
         log_gdp_per_capita = log(gdp_per_capita),
         has_new_cases = ifelse(avg_new_cases > 0, 1, 0),
         log_avg_new_cases = ifelse(has_new_cases, log(avg_new_cases), 0),
         has_new_deaths = ifelse(avg_new_deaths > 0, 1, 0),
         log_avg_new_deaths = ifelse(has_new_deaths, log(avg_new_deaths), 0))

# Create X matrix; center and scale 
X <- model.matrix(~ continent + vaccination_info + has_new_cases +
                    log_population_density + median_age + 
                    log_gdp_per_capita +
                    time*log_avg_new_cases + perc_fully_vaccinated,
                  data = control)
X_scaled <- scale(X[,9:15])
X[,9:15] <- X_scaled
y <- control$int_travel_controls


# Gibb's sampler, uniform priors on cutoffs -------------------------------

# n <- nrow(X)
# J <- length(unique(control$int_travel_controls)) # Number of groups
# ndraws <- 5000
# draws <- matrix(0, nrow = ndraws, ncol = ncol(X) + J-2)
# gamma <- c(-Inf, 0, 0.34, 0.8, 1.44, Inf)
# gamma <- c(-Inf, 0, 1.1, 1.76, 2.7, Inf)
# est_gamma <- 3:(length(gamma) - 1)
# beta <- rep(0.5, ncol(X))
# 
# for (i in 1:ndraws) {
#   z <- rtruncnorm(n, a = gamma[y], b = gamma[y+1], X%*%beta)
#   sd_beta <- solve(t(X)%*%X)
#   beta <- rmvnorm(1, sd_beta %*% t(X)%*%z, sd_beta) %>%
#     t()
#   for (j in est_gamma) {
#     gamma[j] <- runif(1, max(max(z[z >= gamma[j-1] & z < gamma[j]]), gamma[j-1]),
#                       min(min(z[z >= gamma[j] & z < gamma[j+1]]), gamma[j+1]))
#   }
#   draws[i,] <- c(beta, gamma[est_gamma])
# }
# 
# coda::effectiveSize(draws)
# plot(draws[,16], type = 'l')
# plot(draws[,17], type = 'l')
# plot(draws[,18], type = 'l')



# ROUND 2 -----------------------------------------------------------------

doit <- function(ndraws = 10000, warmup = 1000, eps = 0.1) {
  cuts <- c(-Inf, 0, 1, 1.5, 2.5, Inf)
  est.cuts <- 3:(length(cuts)-1)
  nc <- length(est.cuts)
  trans.cuts <- log( c(cuts[est.cuts[1]], diff(cuts[est.cuts])) )
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
    prop.prob <- log(pnorm(upr, mn, 1) - pnorm(lwr, mn, 1))
    sum(prop.prob == -Inf)
    prop.prob[prop.prob == -Inf] <- -1e99
    # Current
    lwr <- cuts[y]
    upr <- cuts[y+1]
    cur.prob <- log(pnorm(upr, mn, 1) - pnorm(lwr, mn, 1))
    sum(cur.prob == -Inf)
    cur.prob[cur.prob == -Inf] <- -1e99
    llike.diff <- sum(prop.prob - cur.prob)
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
  
  cat('Acceptance rate:', accept/ndraws, '\n')
  
  draws <- draws[-(1:warmup),]
  colnames(draws) <- c(colnames(X), 
                       'delta3', 'delta4', 'delta5')
  return(draws)
}

# Run multiple chains
draws <- replicate(4, doit(ndraws = 30000, warmup = 3000, eps = 0.005),
                   simplify = FALSE) %>%
  do.call(rbind, .)

# Do thinning
idx <- seq(4, nrow(draws), by = 4)
draws_init <- draws
draws <- draws[idx,]

# Diagnostic Plots --------------------------------------------------------

Index <- 1:nrow(draws)
p1 <- ggplot(mapping = aes(Index, draws[,1])) +
  geom_line() +
  labs(y = 'Beta_Intercept')
p2 <- ggplot(mapping = aes(Index, draws[,14])) +
  geom_line() +
  labs(y = 'Beta_PropVaccinated')
p3 <- ggplot(mapping = aes(Index, draws[,16])) +
  geom_line() +
  labs(y = 'Delta3')
p4 <- ggplot(mapping = aes(Index, draws[,17])) +
  geom_line() +
  labs(y = 'Delta4')

plot_grid(p1, p2, p3, p4, 
          nrow = 2, byrow = TRUE) %>%
  ggsave('plots/diag-trace.png',
         plot = .,
         device = 'png',
         dpi = 300,
         width = pic_width, 
         height = pic_height,
         units = pic_unit)

ess <- coda::effectiveSize(draws) %>%
  round()

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
                     'Time:Log(Average New Cases)', 'Delta2',
                     'Delta3', 'Delta4')
tbl_ci <- data.frame(coef = tbl_names,
                     est = colMeans(draws) %>%
                       round(3),
                     lower = ci[,1],
                     upper = ci[,2],
                     ess = ess)
rownames(tbl_ci) <- NULL
write_csv(tbl_ci, 'data/tbl_ci.csv')
write.csv(draws, 'data/draws.csv', 
          quote = FALSE, row.names = FALSE)

# Prediction
tw_control <- control %>%
  filter(location == 'Taiwan')
tw_attrs <- tw_control[1,]
X_tw <- data.frame(int = 1,
                     continentAsia = 1,
                     continentEurope = 0,
                     continentNA = 0,
                     continentOceania = 0,
                     continentSA = 0,
                     vaccination_info = c(tw_control$vaccination_info,
                                          rep(1, 6)),
                     has_new_cases = c(tw_control$has_new_cases,
                                       rep(1, 6)),
                     log_population_density = tw_attrs$log_population_density,
                     median_age = tw_attrs$median_age,
                     log_gdp_per_capita = tw_attrs$log_gdp_per_capita,
                     time = c(tw_control$time,
                              23:28),
                     log_avg_new_cases = c(tw_control$log_avg_new_cases,
                                           rep(2, 6)),
                     prop_fully_vaccinated = c(tw_control$perc_fully_vaccinated,
                                               0.4, 0.6, 0.7, 0.75, 0.78, 0.80),
                     time_log_avg_new_cases = c(tw_control$time, 23:28) *
                       c(tw_control$log_avg_new_cases, rep(2, 6)))
contin_mu <- attr(X_scaled, "scaled:center")
contin_sd <- attr(X_scaled, "scaled:scale")
for (j in 9:15) {
  X_tw[,j] <- (X_tw[,j] - contin_mu[j-8]) / contin_sd[j-8]
}
X_tw <- as.matrix(X_tw)

yhat <- matrix(0, nrow = nrow(X_tw), ncol = 3) %>%
  as.data.frame() %>%
  rename(est = V1,
         lower = V2,
         upper = V3) %>%
  mutate(time = 1:nrow(X_tw))
for (i in 1:nrow(X_tw)) {
  prob_open <- apply(draws, 1, function(pars) {
    beta <- pars[1:15]
    # The area under the curve less than this value is the probability
    #  that Taiwan's policies will be at most extreme screening
    delta2 <- pars[16] 
    muhat <- X_tw[i,] %*% beta
    pnorm(delta2, muhat, 1)
  })
  
  yhat[i,1:3] <- c(mean(prob_open), 
                quantile(prob_open, c(0.025, 0.975)))
}

p <- ggplot() +
  geom_polygon(aes(x = c(yhat$time, rev(yhat$time)),
                   y = c(yhat$upper, rev(yhat$lower))),
               fill = 'royalblue', alpha = 0.3) +
  geom_line(aes(time, est),
            data = yhat) +
  labs(x = 'Month (1 = January 2020)',
       y = paste0('P(Y <= 2 | beta, delta)'))

ggsave('plots/results-prediction.png',
       plot = p,
       device = 'png',
       dpi = 300,
       width = pic_width, 
       height = pic_height,
       units = pic_unit)


# Obtain probability draws for April 2022
prob_open <- apply(draws, 1, function(pars) {
  beta <- pars[1:15]
  # The area under the curve less than this value is the probability
  #  that Taiwan's policies will be at most extreme screening
  delta2 <- pars[16] 
  muhat <- X_tw[28,] %*% beta
  pnorm(delta2, muhat, 1)
})


# Plot beta estimates
p <- tbl_ci %>%
  filter( !(coef %in% c('(Intercept)', 'Delta2', 'Delta3', 'Delta4')) ) %>%
  ggplot(aes(x = est, y = coef, xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  geom_errorbar() +
  labs(y = 'Coefficient',
       x = 'Estimate')

ggsave('plots/results-betas.png',
       plot = p,
       device = 'png',
       dpi = 300,
       width = pic_width, 
       height = pic_height,
       units = pic_unit)

