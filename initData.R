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
            start = c(-1.867, 0.5517, 0.3220, 0.6142))
# per a suggestion from https://stats.stackexchange.com/questions/105633/what-to-do-when-a-log-binomial-models-convergence-fails
# replaced 1 in the pmin() function with .9999 to avoid NA/Inf warnings
negLogLik <- function(b) {
  risk <- pmin(.9999, exp(as.matrix(cbind(1, dat[,2:4])) %*% b))
  -sum(dbinom(dat$cancer, 1, risk, log=TRUE))
}
fit <- nlm(negLogLik, p=c(log(mean(dat$cancer)), 0,0,0), hessian=TRUE)
# but careful, it's pretty dependent on starting values
fit <- nlm(negLogLik, p=coef(fit2), hessian=TRUE)
# fiddling with tolerances may help (but not much the way it's written below)
fit <- nlm(negLogLik, p=coef(fit2), hessian=TRUE, gradtol = 1e-12, steptol = 1e-11)
fit <- nlm(negLogLik, p=c(log(mean(dat$cancer)), 0,0,0), hessian=TRUE, gradtol = 1e-12, steptol = 1e-11)
# next steps: can swapping nlm() for optim() help?

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