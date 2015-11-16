# Semi-Parametric-GLARMA
A semi-parametric method for GLARMA type model in Time Series Analysis

The MATLAB function spglarma is used to fit generalized linear autoregressive moving average (GLARMA) models using Pearson residuals whereas the underlying distribution will be estimated from the data with the empirical likelihood approach. The one function therefore will work for discrete, continuous and binary response. 

User can also specify what link they would like to use and currently the options are identity, log, inverse and logit links.



Example Polio 

Cases of Poliomyelitis in the U.S.

Description

This data set gives the monthly number of cases of poliomyelitis in the U.S. for the years 1970–1983 as reported by the Center for Disease Control. The polio data frame has 168 rows and 8 columns.

Format

A csv file containing the following columns:

[:, 1]	Cases	monthly number of cases of poliomyelitis.
[:, 2]	Intcpt	a vector of ones, providing the intercept in the model.
[:, 3]	Trend	a linear trend.
[:, 4]	CosAnnual	cosine harmonics at periods of 12.
[:, 5]	SinAnnual	sine harmonics at periods of 12.
[:, 6]	CosSemiAnnual	cosine harmonics at periods of 6.
[:, 7]	SinSemiAnnual	sine harmonics at periods of 6.

Source

Zeger, S.L (1988) A regression model for time series of counts. Biometrika, 75, 621–629.


% MA(1,2,5), Log link, Pearson Residuals

polio = csvread('./Data/polio.csv',1,2);
Y = polio(:,1);
X = polio(:,2:7);
[delta, maxloglik, fitted, iter, phat, sdhat] = spglarma(Y,X,[],[1,2,5],'log');