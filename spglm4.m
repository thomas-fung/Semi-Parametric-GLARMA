function [beta, maxloglik, fitted, iter, phat, sdhat] = spglm4(y, x, link, tests)
% fits spglm, outputs: 
% estimates of beta
% maximum loglikelihood, 
% number of iterations,
% fitted means, 
% estimated probability mass function phat at each observation
% one row per observation for phat

% the last argument allows you to test hypothesis
% H0: beta[subset] = (some value) via specifying 
% test = X[subset]^T*(some value); set test=0 for joint MLE

n =length(y) ;                      % number of observations
q = size(x,2) ;                     % number of mean parameters

% initial beta value
if  strcmp(link,'id')
    beta0 = [mean(y) zeros(1,q-1)] ;
elseif  strcmp(link,'log')
    beta0 = [log(mean(y)) zeros(1,q-1)] ;
elseif  strcmp(link,'inv')
    beta0 = [1/(mean(y)) zeros(1,q-1)] ;
elseif  strcmp(link,'logit')
    beta0 = [log(mean(y)/(1-mean(y))) zeros(1,q-1)] ;
end;         
p0 = ones(1, n)/n ;                 % initial probabilities
logp0 = log(p0) ;                   % initial log probabilities
b0 = zeros(1, n) ;                  % initial normalizations
theta0 = zeros(1, n) ;              % initial tilts
ell = @(para)loglik4(para, x, y) ;  % passing data to loglik function
conn = @(para)constraints4(para, x, y, link, tests) ; %passing data to constraints
para0 = [beta0 logp0 b0 theta0]' ;  % initial param values

% fit the model
options = optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'TolFun', 1e-8, 'TolCon', 1e-8, 'TolX', 1e-8, 'Algorithm', 'interior-point', 'Display', 'off', 'GradConstr', 'on', 'GradObj', 'on') ;
[param, fvalue, exitflag, output, lambda, grad, hessian]=fmincon(ell, para0, [], [], [], [], [], [],conn, options) ;

% convergence diagnostics
if exitflag == 1||exitflag==2            % convergence indicator
%    fprintf('converged\n');
else
    fprintf('warning: possible boundary issues %d\n', exitflag) ;
end;

% properties of the fitted model
beta = param(1:q) ;                 % fitted beta
maxloglik = -fvalue ;               % loglikeliood at maximim
iter = output.iterations ;          % number of iterations
fncount = output.funcCount ;

% fitted mean values
if  strcmp(link,'id')
    fitted = x*beta + tests; 
elseif  strcmp(link,'log')
    fitted = exp(x*beta + tests);
elseif  strcmp(link,'inv')
    fitted = 1./(x*beta + tests);
elseif  strcmp(link,'logit')
    fitted = exp(x*beta+tests)./(1+exp(x*beta+tests));
end;                 

% estimated distribution and standard deviation at each observation
% if requested
if nargout > 4
    p = exp(param(q+1:q+n)) ;           % fitted p
    b = param(q+n+1:q+2*n) ;            % fitted b
    theta = param(q+2*n+1: q+3*n) ;     % fitted theta
    phat = ones(n, n);
    Varhat = ones(n,1);
    for i=1:n
        phat(i,:) = p.*exp(b(i)+theta(i)*y);
        Varhat(i) = phat(i,:)*((y-fitted(i)).^2) ;
    end;
    sdhat = sqrt(Varhat) ;
end;

end
