function [l, gradl] = loglik4(para, x, y) 
% computes the negative loglikelihood
% y is a vertical n-vector, x is a n by p matrix

n = length(y) ;
q = size(x,2) ;
%beta = para(1:q) ;             %mean parameters
logp = para(q+1:q+n) ;          %log probability vector
b = para(q+n+1: q+2*n) ;        %normalizing constants
theta =para(q+2*n+1:q+3*n) ;    %tilt values
l = -sum(logp + b + theta.*y) ;
if nargout > 1 
    gradl = [zeros(q,1); -ones(n,1) ; -ones(n,1) ; -y]; 
    % vertical vector of gradients
end;
end 