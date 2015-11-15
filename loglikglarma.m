function [l, gradl] = loglikglarma(para, Y, X, phi_lags, theta_lags)

% ver 1.0.0
% used in conjunction with spglarma.m ver 1.0.0

% This function acts as the objective function of the sp glarma algorithm, outputs:
% l = negative loglikelihood
% gradl = gradients of the objective function with respect to the parameters

% usage [GL] = loglikglarma(Y, X, phi_lags, theta_lags,link)
% para = Parameters [beta' phitheta Z e' logp b xi]
% Y = Response Variable (e.g. Non-negative Integers) n by 1
% X = Matrix of covariates n by p
% phi_lags = AR X-th Order (1,2,4,6...) in a row vector; e.g. [1 2 4]
% theta_lags = MA X-th Order (1,7...) in a row vector; e.g. [1 2 5]
% link = 'id', 'log', 'inv' or 'logit'

% computes the negative loglikelihood
% y is a vertical n-vector, x is a n by p matrix
n = length(Y) ;
rsq = size(X,2)+ length(phi_lags)+ length(theta_lags); 
logp = para(rsq+2*n+1:rsq+3*n) ;          %log probability vector
b = para(rsq+3*n+1: rsq+4*n) ;        %normalizing constants
xi =para(rsq+4*n+1:rsq+5*n) ;    %tilt values
l = -sum(logp + b + xi.*Y') ;
if nargout > 1 
    gradl = [zeros(rsq+2*n,1); -ones(n,1) ; -ones(n,1) ; -Y]; 
    % vertical vector of gradients
end;
end 