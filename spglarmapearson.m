function [delta, maxloglik, fitted, iter, phat, sdhat] = spglarmapearson(Y, X, phi_lags, theta_lags,link,tests)

n = length(Y); % number of observation
s = length(phi_lags);
q = length(theta_lags);
r = size(X,2);
rsq = r+s+q;

% [beta0, ~, fitted, ~, ~, sdhat] = spglm4(Y,X,link, 0); %initial beta estimates
if  strcmp(link,'id')
    beta0 = [mean(Y); zeros(r-1,1)] ;
elseif  strcmp(link,'log')
    beta0 = [log(mean(Y)); zeros(r-1,1)] ;
elseif  strcmp(link,'inv')
    beta0 = [1/(mean(Y)); zeros(r-1,1)] ;
elseif  strcmp(link,'logit')
    beta0 = [log(mean(Y)/(1-mean(Y))); zeros(r-1,1)] ;
end
phitheta0 = zeros(1,s+q);  % initial phi and theta estimates
Z0 = zeros(1, n); %initial Z estimates;
p0 = ones(1, n)/n ;                 % initial probabilities
logp0 = log(p0) ;                   % initial log probabilities
b0 = zeros(1, n) ;                  % initial normalizations
xi0 = zeros(1, n) ;              % initial tilts

ell = @(para)loglikglarma(para, Y, X, phi_lags, theta_lags) ;  % passing data to loglik function
conn = @(para)constraintsglarmapearson(para, Y, X, phi_lags, theta_lags,link) ; %passing data to constraints
para0 = [beta0' phitheta0 Z0 logp0 b0 xi0]; % initial param values
    

% fit the model
%options = optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'TolFun', 1e-8, 'TolCon', 1e-8, 'TolX', 1e-8, 'Algorithm', 'interior-point', 'Display', 'off', 'GradConstr', 'on', 'GradObj', 'on',) ;
options = optimoptions('fmincon','MaxFunEvals',1e5, 'MaxIter', 1e5, 'TolFun', 1e-8, ...
    'TolCon', 1e-8, 'Algorithm', 'sqp', 'StepTolerance', 1e-8, ...
    'Display', 'off', 'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true,...
    'UseParallel',false);

if isempty(tests)
    [param, fvalue, exitflag, output, ~, ~, ~]=fmincon(ell, para0, [], [], [], [], [], [], conn, options) ;
else
    lbrsq = -inf*ones(1,rsq);
    uprsq = inf*ones(1,rsq);
    lbrsq(~isnan(tests)) = tests(~isnan(tests));
    uprsq(~isnan(tests)) = tests(~isnan(tests));
    lb = [lbrsq, -inf*ones(1,4*n)];
    ub = [uprsq, inf*ones(1,4*n)];
    [param, fvalue, exitflag, output, ~, ~, ~]=fmincon(ell, para0, [], [], [], [], lb, ub, conn, options);
end

%options = optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'TolFun', 1e-8, 'TolCon', 1e-8, 'TolX', 1e-8, 'Algorithm', 'interior-point', 'Display', 'off', 'GradConstr', 'on', 'GradObj', 'on','UseParallel','always');
%[param, fvalue, exitflag, output, ~, ~, ~]=fmincon(ell, para0, [], [], [], [], [], [], conn, options) ;
%[param, fvalue, exitflag, output, lambda, grad, hessian]=fmincon(ell, para0, [], [], [], [], [], [], conn, options) ;

% convergence diagnostics
if exitflag == 1||exitflag==2            % convergence indicator
%    fprintf('converged\n');
else
    fprintf('warning: possible boundary issues %d\n', exitflag) ;
end

% properties of the fitted model


maxloglik = -fvalue ;               % loglikeliood at maximim
iter = output.iterations ;          % number of iterations
fncount = output.funcCount ;


delta = param(1:rsq)';
beta = param(1:r)';
Z = param(rsq+1:rsq+n)';
% fitted mean values
if  strcmp(link,'id')
    fitted = X*beta+ Z; 
elseif  strcmp(link,'log')
    fitted = exp(X*beta+ Z);
elseif  strcmp(link,'inv')
    fitted = 1./(X*beta +Z);
elseif  strcmp(link,'logit')
    fitted = exp(X*beta+Z)./(1+exp(X*beta+Z));
end                 

% estimated distribution and standard deviation at each observation
% if requested
if nargout > 4
    p = exp(param(rsq+n+1:rsq+2*n)) ;           % fitted p
    b = param(rsq+2*n+1:rsq+3*n) ;            % fitted b
    xi = param(rsq+3*n+1: rsq+4*n) ;     % fitted theta
    phat = ones(n, n);
    Varhat = ones(n,1);
    for i=1:n
        phat(i,:) = p.*exp(b(i)+xi(i)*Y');
        Varhat(i) = phat(i,:)*((Y-fitted(i)).^2) ;
    end
    sdhat = sqrt(Varhat) ;
end

end
