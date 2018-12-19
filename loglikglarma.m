function [l, gradl] = loglikglarma(para, Y, X, phi_lags, theta_lags) 
% computes the negative loglikelihood
% y is a vertical n-vector, x is a n by p matrix
n = length(Y) ;
rsq = size(X,2)+ length(phi_lags)+ length(theta_lags); 
logp = para(rsq+n+1:rsq+2*n) ;          %log probability vector
b = para(rsq+2*n+1: rsq+3*n) ;        %normalizing constants
xi =para(rsq+3*n+1:rsq+4*n) ;    %tilt values
l = -sum(logp + b + xi.*Y') ;
if nargout > 1 
    gradl = [zeros(rsq+n,1); -ones(n,1) ; -ones(n,1) ; -Y]; 
    % vertical vector of gradients
end
end 