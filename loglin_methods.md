
``` r
library(tidyverse)
```

``` r
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
      prob = 1 / (1 + exp(2 - .$PSA - .$famhx - .$race))
    )
  )
```

### Log binomial fit

``` r
glm(cancer ~ PSA + famhx + race, data = dat, family = binomial(link = "log"))
```

    ## Error: no valid set of coefficients has been found: please supply starting values

### Poisson fit

``` r
fit2 <- glm(cancer ~ PSA + famhx + race, data = dat, family = poisson)
fit2
```

    ## 
    ## Call:  glm(formula = cancer ~ PSA + famhx + race, family = poisson, 
    ##     data = dat)
    ## 
    ## Coefficients:
    ## (Intercept)          PSA        famhx         race  
    ##     -1.8671       0.5517       0.3220       0.6142  
    ## 
    ## Degrees of Freedom: 999 Total (i.e. Null);  996 Residual
    ## Null Deviance:       719.4 
    ## Residual Deviance: 650   AIC: 1244

### Poisson fit with lme4

``` r
fit3 <- lme4::glmer(cancer ~ PSA + famhx + race + (1 | id),
  data = dat, family = poisson
)
fit3
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: cancer ~ PSA + famhx + race + (1 | id)
    ##    Data: dat
    ##       AIC       BIC    logLik  deviance  df.resid 
    ## 1246.0458 1270.5846 -618.0229 1236.0458       995 
    ## Random effects:
    ##  Groups Name        Std.Dev. 
    ##  id     (Intercept) 1.579e-08
    ## Number of obs: 1000, groups:  id, 1000
    ## Fixed Effects:
    ## (Intercept)          PSA        famhx         race  
    ##     -1.8671       0.5517       0.3220       0.6142  
    ## convergence code 0; 1 optimizer warnings; 0 lme4 warnings

### Log binomial fit with starting values from Poisson

In theory, this should work with `coef(fit2)` in the start param, but
that throws an error because these coefficients are either outside of or
too close to the feasible region boundary where coefficients are bound
to `[0, 1]`.

``` r
coef(fit2)
```

    ## (Intercept)         PSA       famhx        race 
    ##  -1.8670565   0.5517099   0.3220246   0.6142025

``` r
fit4 <- glm(cancer ~ PSA + famhx + race,
  data = dat,
  family = binomial(link = "log"),
  start = coef(fit2)
)
```

    ## Error: cannot find valid starting values: please specify some

## negLogLik

Per a suggestion from
[stats.stackexchange](https://stats.stackexchange.com/questions/105633/what-to-do-when-a-log-binomial-models-convergence-fails),
but replacing 1 in the `pmin()` function with .9999 to avoid `NA`/`Inf`
warnings.

``` r
# https://stats.stackexchange.com/a/321407/20236
negLogLik <- function(b) {
  risk <- pmin(.9999, exp(as.matrix(cbind(1, dat[, 2:4])) %*% b))
  -sum(dbinom(dat$cancer, 1, risk, log = TRUE))
}
(fit <- nlm(negLogLik, p = c(log(mean(dat$cancer)), 0, 0, 0), hessian = TRUE))
```

    ## $minimum
    ## [1] 553.7325
    ## 
    ## $estimate
    ## [1] -1.7942923  0.4648994  0.3181380  0.5758614
    ## 
    ## $gradient
    ## [1]  0.0001501638 -0.0001085709  0.0000172804  0.0000163709
    ## 
    ## $hessian
    ##           [,1]       [,2]     [,3]      [,4]
    ## [1,] 600.63148  724.04328 35.09289 419.11532
    ## [2,] 724.04328 1400.77481 22.26708 592.21486
    ## [3,]  35.09289   22.26708 35.09289  25.17237
    ## [4,] 419.11532  592.21486 25.17237 419.11532
    ## 
    ## $code
    ## [1] 1
    ## 
    ## $iterations
    ## [1] 14

but careful, it’s pretty dependent on starting values

``` r
(fit <- nlm(negLogLik, p = coef(fit2), hessian = TRUE))
```

    ## $minimum
    ## [1] 565.6537
    ## 
    ## $estimate
    ## [1] -1.8741224  0.5567255  0.2605088  0.6422139
    ## 
    ## $gradient
    ## [1]  2.729762e-05  3.637979e-05 -1.227818e-05  4.786216e-05
    ## 
    ## $hessian
    ##           [,1]       [,2]     [,3]      [,4]
    ## [1,] 542.25621  565.12544 33.62811 325.25674
    ## [2,] 565.12544 1013.58905 22.48684 303.14848
    ## [3,]  33.62811   22.48684 33.62811  25.01379
    ## [4,] 325.25674  303.14848 25.01379 325.25672
    ## 
    ## $code
    ## [1] 1
    ## 
    ## $iterations
    ## [1] 18

fiddling with tolerances may help (but not much the way it’s written
below)

``` r
(fit <- nlm(negLogLik, p = coef(fit2), hessian = TRUE, gradtol = 1e-12, steptol = 1e-11))
```

    ## $minimum
    ## [1] 565.6537
    ## 
    ## $estimate
    ## [1] -1.8741205  0.5567252  0.2605089  0.6422125
    ## 
    ## $gradient
    ## [1] 0.000000e+00 0.000000e+00 5.684342e-10 0.000000e+00
    ## 
    ## $hessian
    ##           [,1]       [,2]     [,3]      [,4]
    ## [1,] 542.25779  565.12843 33.62819 325.25710
    ## [2,] 565.12843 1013.59710 22.48689 303.14882
    ## [3,]  33.62819   22.48689 33.62818  25.01384
    ## [4,] 325.25710  303.14882 25.01384 325.25710
    ## 
    ## $code
    ## [1] 2
    ## 
    ## $iterations
    ## [1] 38

``` r
(fit <- nlm(negLogLik, p = c(log(mean(dat$cancer)), 0, 0, 0), hessian = TRUE, gradtol = 1e-12, steptol = 1e-11))
```

    ## $minimum
    ## [1] 553.7325
    ## 
    ## $estimate
    ## [1] -1.7942930  0.4648998  0.3181379  0.5758616
    ## 
    ## $gradient
    ## [1]  5.810134e-05 -9.208634e-06  7.958079e-07  1.023182e-06
    ## 
    ## $hessian
    ##           [,1]       [,2]     [,3]      [,4]
    ## [1,] 600.63242  724.04612 35.09284 419.11635
    ## [2,] 724.04612 1400.78190 22.26706 592.21767
    ## [3,]  35.09284   22.26706 35.09284  25.17233
    ## [4,] 419.11635  592.21767 25.17233 419.11635
    ## 
    ## $code
    ## [1] 2
    ## 
    ## $iterations
    ## [1] 24

next steps: can swapping nlm() for optim() help?

## Logbin

``` r
library(logbin)
(fit <- logbin(cancer ~ PSA + famhx + race, data = dat))
```

    ## 
    ## Call:  logbin(formula = cancer ~ PSA + famhx + race, data = dat)
    ## 
    ## Coefficients:
    ## (Intercept)          PSA        famhx         race  
    ##  -1.717e+00    3.677e-01    7.012e-07    5.767e-01  
    ## 
    ## Degrees of Freedom: 999 Total (i.e. Null);  996 Residual
    ## Null Deviance:       1210 
    ## Residual Deviance: 1115  AIC: 1123

It doesn’t compute cov matrix to enable CI estimation\!

We could use inverted test-based limits [according to
Agresti](http://statmath.wu.ac.at/research/talks/resources/slidesagresti_confidence.pdf)
and then simulate coverage review of test-based CIs.

## Example from Spiegelman and Hertzmark AJE 2005

``` r
dat <- tibble(
  id = 1:192,
  death = c(rep(1, 54), rep(0, 138)),
  stage = c(
    rep("Stage I", 7), rep("Stage II", 26), rep("Stage III", 21),
    rep("Stage I", 60), rep("Stage II", 70), rep("Stage III", 8)
  ),
  receptor = c(
    rep("Low", 2), rep("High", 5), rep("Low", 9), rep("High", 17),
    rep("Low", 12), rep("High", 9), rep("Low", 10), rep("High", 50),
    rep("Low", 13), rep("High", 57), rep("Low", 2), rep("High", 6)
  )
)
dat
```

    ## # A tibble: 192 x 4
    ##       id death stage    receptor
    ##    <int> <dbl> <chr>    <chr>   
    ##  1     1     1 Stage I  Low     
    ##  2     2     1 Stage I  Low     
    ##  3     3     1 Stage I  High    
    ##  4     4     1 Stage I  High    
    ##  5     5     1 Stage I  High    
    ##  6     6     1 Stage I  High    
    ##  7     7     1 Stage I  High    
    ##  8     8     1 Stage II Low     
    ##  9     9     1 Stage II Low     
    ## 10    10     1 Stage II Low     
    ## # ... with 182 more rows

log binomial fails to
converge

``` r
glm(death ~ stage + receptor, data = dat, family = binomial(link = "log"))
```

    ## Error: no valid set of coefficients has been found: please supply starting values

Poisson converges

``` r
(fit <- glm(death ~ stage + receptor, data = dat, family = poisson))
```

    ## 
    ## Call:  glm(formula = death ~ stage + receptor, family = poisson, data = dat)
    ## 
    ## Coefficients:
    ##    (Intercept)   stageStage II  stageStage III     receptorLow  
    ##        -2.3658          0.9246          1.7772          0.4891  
    ## 
    ## Degrees of Freedom: 191 Total (i.e. Null);  188 Residual
    ## Null Deviance:       137 
    ## Residual Deviance: 110.3     AIC: 226.3

``` r
exp(coef(fit)) #' estimates match https://academic.oup.com/aje/article/162/3/199/171116
```

    ##    (Intercept)  stageStage II stageStage III    receptorLow 
    ##     0.09387241     2.52074190     5.91337206     1.63077519

lme4
Poisson

``` r
(fit3 <- lme4::glmer(death ~ as.factor(stage) + receptor + (1 | id), data = dat, family = poisson))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: death ~ as.factor(stage) + receptor + (1 | id)
    ##    Data: dat
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  228.2920  244.5795 -109.1460  218.2920       187 
    ## Random effects:
    ##  Groups Name        Std.Dev.
    ##  id     (Intercept) 0       
    ## Number of obs: 192, groups:  id, 192
    ## Fixed Effects:
    ##               (Intercept)   as.factor(stage)Stage II  
    ##                   -2.3658                     0.9246  
    ## as.factor(stage)Stage III                receptorLow  
    ##                    1.7772                     0.4891  
    ## convergence code 0; 1 optimizer warnings; 0 lme4 warnings

using Poisson as starting
values

``` r
(fit2 <- glm(death ~ as.factor(stage) + receptor, data = dat, family = binomial(link = "log"), start = c(-2.3, .92, 1.78, .489)))
```

    ## 
    ## Call:  glm(formula = death ~ as.factor(stage) + receptor, family = binomial(link = "log"), 
    ##     data = dat, start = c(-2.3, 0.92, 1.78, 0.489))
    ## 
    ## Coefficients:
    ##               (Intercept)   as.factor(stage)Stage II  
    ##                   -2.3521                     0.9314  
    ## as.factor(stage)Stage III                receptorLow  
    ##                    1.7695                     0.4436  
    ## 
    ## Degrees of Freedom: 191 Total (i.e. Null);  188 Residual
    ## Null Deviance:       228.1 
    ## Residual Deviance: 185.9     AIC: 193.9

``` r
exp(coef(fit2))
```

    ##               (Intercept)  as.factor(stage)Stage II 
    ##                0.09516633                2.53815743 
    ## as.factor(stage)Stage III               receptorLow 
    ##                5.86805395                1.55832434
