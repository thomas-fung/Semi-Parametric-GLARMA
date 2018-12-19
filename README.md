# Semi-Parametric-GLARMA
## A semi-parametric method for GLARMA type model in Time Series Analysis (Matlab Version)

The MATLAB function `spglarma` is used to fit generalized linear autoregressive moving average (GLARMA) models using Pearson or Score-type residuals whereas the underlying distribution will be estimated from the data with the empirical likelihood approach of [Huang and Fung (2016)](https://arxiv.org/abs/1603.02802). 

Currently only the 'log' link option is fully tested. There are plans to implement and test the identity, inverse and logit links options.

An `R` package `spglarma` is currently under development and will be released on github and cran shortly to accompany an updated version of the paper.

## Installation
This implmentation is based around using the `fmincon` function of `Matlab`'s optimization toolbox. Saved the following m-files into a working directory: `spglarma`, `spglarmapearson`, `spglarmascore`, `loglikglarma`, `constraintsglarmapearson`, `constraintsglarmascore`, `plotspglarma`, `SPPIT_hist`, `spglm4`, `loglik4`, `constraints4` and then you are good to go.  

## Citation

If you use these code to analyse your data, please use the following citation:

- Huang, A. and Fung, T. (2016). Semiparametric generalized linear models for time-series data. arXiv:1603.02802.
 
## Polio Example 

Cases of Poliomyelitis in the U.S.

### Description

This data set gives the monthly number of cases of poliomyelitis in the U.S. for the years 1970–1983 as reported by the Center for Disease Control. The polio data frame has 168 rows and 8 columns.

### Format

A csv file containing the following columns:

- [:, 1]	Cases	monthly number of cases of poliomyelitis.
- [:, 2]	Intcpt,	a vector of ones, providing the intercept in the model.
- [:, 3]	Trend	a linear trend. 
- [:, 4]	CosAnnual	cosine harmonics at periods of 12.
- [:, 5]	SinAnnual	sine harmonics at periods of 12.
- [:, 6]	CosSemiAnnual	cosine harmonics at periods of 6.
- [:, 7]	SinSemiAnnual	sine harmonics at periods of 6.

### Source

Zeger, S.L (1988) A regression model for time series of counts. Biometrika, 75, 621–629.

### Sample code
To fit a model with MA(1,2,5), Log link, Pearson Residuals
```s
polio = csvread('./Data/polio.csv',1,2);
Y = polio(:,1);
X = polio(:,2:7);
[delta, maxloglik, fitted, iter, phat, sdhat] = spglarma(Y,X,[],[1,2,5],'log');
```

More details of the example can be found in `GLARMAExample_Polio_github.m`.
