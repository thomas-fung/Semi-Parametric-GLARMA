function [delta, maxloglik, fitted, iter, phat, sdhat] = spglarma(Y, X, phi_lags, theta_lags,link,tests,residualtype)

% Main file for fitting the semi-parametric glarma model 
% with the fmincon and sqp algorithm 

%usage [GL] = spglarma(Y, X, phi_lags, theta_lags,link)
% Y = Response Variable (e.g. Non-negative Integers) n by 1
% X = Matrix of covariates n by p 
% phi_lags = AR X-th Order (1,2,4,6...) in a row vector 
% theta_lags = MA X-th Order (1,7...) in a row vector
% link = 'log' (need to check for 'id', 'inv' or 'logit' implementation)

% the last argument allows you to test hypothesis
% H0: beta[subset] = (some value) via specifying 
% test = X[subset]^T*(some value); set test=0 for joint MLE

% adding a constant column for the intercept if it is not part of X already
% if not(all(X(:,1)==1))
%     X=[ones(size(X,1),1),X];
% end


% the last argument allows you to test hypothesis
% H0: beta[subset] = (some value) via specifying 
% test = X[subset]^T*(some value); set test=0 for joint MLE

% 'Pearson' for pearson residual
% 'Score' for Score-type residual
switch residualtype 
    case 'Pearson' 
        [delta, maxloglik, fitted, iter, phat, sdhat] = spglarmapearson(Y, X, phi_lags, theta_lags,link,tests);
    case 'Score'
        [delta, maxloglik, fitted, iter, phat, sdhat] = spglarmascore(Y, X, phi_lags, theta_lags,link,tests);

    otherwise
        warning('Unexpected residual type')
end


end
