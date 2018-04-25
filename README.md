jagsRun
====

[![Build Status](https://travis-ci.org/caseyyoungflesh/jagsRun.svg?branch=master)](https://travis-ci.org/caseyyoungflesh/jagsRun)

`jagsRun` is an R package used run JAGS in parallel and produce summary output.

The package contains one function:

- `jagsRun` - run JAGS in parallel


Installation
------------

You can install the  development version from Github with:
```{r}
install.packages('devtools')
devtools::install_github('caseyyoungflesh/jagsRun')
```

If there are firewall issues or issues with libcurl when using devtools the following can be run from bash terminal to install package in specified location (with cluster use):
```
wget --no-check-certificate https://github.com/caseyyoungflesh/jagsRun/archive/master.tar.gz
R CMD INSTALL master.tar.gz
```
