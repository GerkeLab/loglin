
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Easy R Code for Convergent Log-Linear Models

<!--can probably cite https://www.ncbi.nlm.nih.gov/pubmed/3509965 in the first sentence but need to read it to confirm-->

## Introduction

Reporting risk ratios instead of odds ratios is often preferable for
interpretability. In a regression framework, risk ratios are most
efficiently estimated through log-binomial generalized linear models;
unfortunately, statistical software algorithms often fail to converge in
the log-binomial setting (1). A previous publication in this journal
introduced a SAS macro that can, in some settings, leverage parameter
estimates from a modified Poisson regression to guide the log-binomial
fitting algorithm towards convergence (2). When this approach fails to
fit a log-binomial model, parameter estimates from the modified Poisson
regression are reported, though such estimates are known to exhibit
reduced statistical efficiency.

Here, we present …

## Example

Example data, procedure, etc.

## References

<div id="refs" class="references">

<div id="ref-Williamson2013-Logbinomial">

1\. Williamson T, Eliasziw M, Fick GH. Log-binomial models: Exploring
failed convergence. *Emerg Themes Epidemiol*. 2013;10:14. 

</div>

<div id="ref-Spiegelman2005-Easy">

2\. Spiegelman D, Hertzmark E. Easy SAS Calculations for Risk or
Prevalence Ratios and Differences. *Am J Epidemiol*.
2005;162(3):199–200. 

</div>

</div>
