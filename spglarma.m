function [delta, maxloglik, fitted, iter, phat, sdhat] = spglarma(Y, X, phi_lags, theta_lags,link)

% version 1.0.0

% fits sp glarma with pearson residuals, outputs:
% delta = estimates of beta, ARMA coefficients
% maxloglik = maximum loglikelihood,
% fitted = fitted values,
% iter = number of iterations till convergence,
% phat = estimated probability mass function phat at each observation; one row per observation for phat
% sdhat = the estimated standard deviation

% usage [GL] = GLARMA(Y, X, phi_lags, theta_lags,link)
% Y = Response Variable (e.g. Non-negative Integers) n by 1
% X = Matrix of covariates n by p 
% phi_lags = AR X-th Order (1,2,4,6...) in a row vector; e.g. [1 2 4]
% theta_lags = MA X-th Order (1,7...) in a row vector; e.g. [1 2 5]
% link = 'id', 'log', 'inv' or 'logit'

% adding a constant column for the intercept if it is not part of X already
if not(all(X(:,1)==1))
    X=[ones(size(X,1),1),X];
end

% getting the size of various vectors/matrices
n = length(Y); % number of observation
s = length(phi_lags);
q = length(theta_lags);
r = size(X,2);
rsq = r+s+q;


% initial estimates
[beta0, ~, fitted, ~, ~, sdhat] = spglm4(Y,X,link, 0); % initial beta estimates provided by AH's spglm algorithm
e0 = (Y - fitted)./sdhat; % initial pearson residuals estimates using beta only
phitheta0 = zeros(1,s+q);  % initial phi and theta estimates
Z0 = zeros(1, n); %initial Z estimates;
p0 = ones(1, n)/n ;                 % initial probabilities
logp0 = log(p0) ;                   % initial log probabilities
b0 = zeros(1, n) ;                  % initial normalizations
xi0 = zeros(1, n) ;              % initial tilts

% initiates the objective & constraint functions
ell = @(para)loglikglarma(para, Y, X, phi_lags, theta_lags) ;  % passing data to loglik function
conn = @(para)constraintsglarma(para, Y, X, phi_lags, theta_lags,link) ; %passing data to constraints
para0 = [beta0' phitheta0 Z0 e0' logp0 b0 xi0]; % initial param values
    

% fit the model
options = optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'TolFun', 1e-8, 'TolCon', 1e-8, 'TolX', 1e-8, 'Algorithm', 'interior-point', 'Display', 'off', 'GradConstr', 'off', 'GradObj', 'on') ;
[param, fvalue, exitflag, output, ~, ~, ~]=fmincon(ell, para0, [], [], [], [], [], [], conn, options) ;

% convergence diagnostics
if exitflag == 1||exitflag==2            % convergence indicator
%    fprintf('converged\n');
else
    fprintf('warning: possible boundary issues %d\n', exitflag) ;
end;

% properties of the fitted model
maxloglik = -fvalue ;               % loglikeliood at maximim
iter = output.iterations ;          % number of iterations
fncount = output.funcCount ;


delta = param(1:rsq)';
beta = param(1:r)';
Z = param(rsq+1:rsq+n)';
% fitted mean values
if  strcmp(link,'id')
    fitted = X*beta + Z; 
elseif  strcmp(link,'log')
    fitted = exp(X*beta + Z);
elseif  strcmp(link,'inv')
    fitted = 1./(X*beta +Z);
elseif  strcmp(link,'logit')
    fitted = exp(X*beta+Z)./(1+exp(X*beta+Z));
end;                 

% estimated distribution and standard deviation at each observation
% if requested
if nargout > 4
    p = exp(param(rsq+2*n+1:rsq+3*n)) ;           % fitted p
    b = param(rsq+3*n+1:rsq+4*n) ;            % fitted b
    xi = param(rsq+4*n+1: rsq+5*n) ;     % fitted theta
    phat = ones(n, n);
    Varhat = ones(n,1);
    for i=1:n
        phat(i,:) = p.*exp(b(i)+xi(i)*Y');
        Varhat(i) = phat(i,:)*((Y-fitted(i)).^2) ;
    end;
    sdhat = sqrt(Varhat) ;
end;

end
