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

#### Example from Spiegelman and Hertzmark AJE 2005
dat <- tibble(id = 1:192,
              death = c(rep(1, 54), rep(0, 138)),
              stage = c(rep("Stage I", 7), rep("Stage II", 26), rep("Stage III", 21),
                        rep("Stage I", 60), rep("Stage II", 70), rep("Stage III", 8)),
              receptor = c(rep("Low", 2), rep("High", 5), rep("Low", 9), rep("High", 17),
                           rep("Low", 12), rep("High", 9), rep("Low", 10), rep("High", 50),
                           rep("Low", 13), rep("High", 57), rep("Low", 2), rep("High", 6)))

# log binomial fails to converge
glm(death ~ stage + receptor, data = dat, family = binomial(link = "log"))
# Poisson converges
fit <- glm(death ~ stage + receptor, data = dat, family = poisson)
exp(coef(fit)) # estimates match https://academic.oup.com/aje/article/162/3/199/171116
# lme4 Poisson
fit3 <- lme4::glmer(death ~ as.factor(stage) + receptor + (1 | id), data = dat, family=poisson)
# using Poisson as starting values
fit2 <- glm(death ~ as.factor(stage) + receptor, data = dat, family = binomial(link = "log"), start = c(-2.3, .92, 1.78, .489))
exp(coef(fit2))