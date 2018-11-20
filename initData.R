library(tidyverse)

n <- 1000
set.seed(8675309)
dat <- tibble(
  id = 1:n,
  PSA = rexp(n, 2),
  famhx = rbinom(n, size = 1, prob = .05),
  race = rbinom(n, size = 1, prob = .4)
) %>%
  add_column(
    cancer = rbinom(
      n,
      size = 1,
      prob = 1/(1+exp(2 - .$PSA - .$famhx - .$race))
    )
  )

# log binomial fit
glm(cancer ~ PSA + famhx + race, data = dat, family = binomial(link = "log"))
# Poisson fit
fit2 <- glm(cancer ~ PSA + famhx + race, data = dat, family = poisson)
# Poisson fit with lme4
fit3 <- lme4::glmer(cancer ~ PSA + famhx + race + (1 | id),
                    data = dat, family = poisson)
# log binomial fit with starting values from Poisson
# this should work with coef(fit2) in the start param, but that throws errors
fit4 <- glm(cancer ~ PSA + famhx + race, data = dat,
            family = binomial(link = "log"),
            start = c(-1.87, 0.55, 0.322, 0.614))
