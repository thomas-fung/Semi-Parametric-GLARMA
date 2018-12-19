function [c, ceq, gradc, gradceq] = constraintsglarmascore(para, Y, X, phi_lags, theta_lags, link)
% computes the mean and normalization constraints
% y is a vertical n-vector, x is a n by p matrix

n = length(Y) ;
r = size(X,2) ;
s = length(phi_lags); % number of AR component 
q = length(theta_lags); % number of MA component
beta = para(1:r)' ; %mean parameters
rsq = r+s+q;
if s>0 && q>0 
    phi = para(r+1:r+s); %AR parameters
    theta = para(r+s+1:rsq); % MA parameters
elseif s>0 && q<=0
    phi = para(r+1:r+s); % AR parameters/ no MA parameters
else
    theta = para(r+1:r+q); % MA parameters/ no AR parameters
end
Z = para(rsq+1:rsq+n)'; % ARMA component parameters
p = exp(para(rsq+n+1:rsq+2*n)) ;        %probability vector
b = para(rsq+2*n+1:rsq+3*n) ;        %normalizing constants
xi =para(rsq+3*n+1:rsq+4*n) ;    %tilt values
if s>0 && q <=0 
    msq = phi_lags(s);
elseif s<=0 && q>0
    msq = theta_lags(q);
else
    msq = max(phi_lags(s),theta_lags(q));
end



% mean values under beta and test constraints
if  strcmp(link,'id')
    mu = X*beta+ Z; 
elseif  strcmp(link,'log')
    mu = exp(X*beta + Z);
elseif  strcmp(link,'inv')
    mu = 1./(X*beta+ Z);
elseif  strcmp(link,'logit')
    mu = exp(X*beta+Z)./(1+exp(X*beta+Z)) ; 
end

% to store mean and norm constraints
% can be written more efficiently without for loops
mmu = zeros(1,n);
nnorm = zeros(1,n) ;
phat = ones(n, n);
Varhat = ones(n,1);
    for i=1:n
        phat(i,:) = p.*exp(b(i)+xi(i)*Y');
        Varhat(i) = phat(i,:)*((Y-mu(i)).^2) ;
    end

e = (Y-mu)./Varhat; % score residuals

nmsq = n+msq;
%Setting up the ARMA constraints;
Ztime = [zeros(msq,1); Z];
Zc = zeros(n,1);
etime = [zeros(msq,1);e];

if(s > 0) 
              for i = 1:s
                    Zc = Zc + phi(i) .* (Ztime(msq+1- phi_lags(i):nmsq-phi_lags(i))+etime(msq+1-phi_lags(i):nmsq-phi_lags(i)));
                end
            end
       
            if(q > 0) 
                for i = 1:q
            		Zc = Zc + theta(i) * etime(msq+1-theta_lags(i):nmsq - theta_lags(i));
                end
            end


% compute mean, norm constraints and their gradients (if requested)
if nargout > 2
    gradoutm = zeros(rsq+5*n, n) ; % column i = mean constraint i
    gradoutn = zeros(rsq+5*n, n) ; % column i = norm constraint i
    gradouta = zeros(rsq+5*n, n) ; % column i = arma constraint i
    
    for i=1:n
        % precompute each of these quantities
        pexponenti = p.*exp(b(i)+xi(i)*Y') ;
        ypexponenti = Y'.*p.*exp(b(i)+xi(i)*Y') ;
        y2pexponenti = Y'.^2.*p.*exp(b(i)+xi(i)*Y') ;
        
        % mean and norm constraints
        mmu(i)= sum(ypexponenti) ;
        nnorm(i)= sum(pexponenti) ;
        
        % compute gradients of mean constraint
        if  strcmp(link,'id')
            gradoutm(:,i) = [X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'log')
            gradoutm(:,i) = [mu(i)*X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'inv')
            gradoutm(:,i) = [-mu(i)^2*X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'logit')
            gradoutm(:,i) = [mu(i)*(1-mu(i))*X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        end;
        
        % compute gradient of norm constraint
        gradoutn(:,i) = [zeros(q,1) ; -pexponenti ; [zeros(1,i-1) -nnorm(i) zeros(1, n-i)]' ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]'] ;

	% compute gradient of the arma constraint

	% compute graduent of the residual constraint
	
    end
    
    outm = transpose(mu)-mmu ; % outm is a horizontal vector
    outn = 1 - nnorm ;      % outn is a horiztonal vector
    outarma = transpose(Z)- transpose(Zc);            % arma constraints
    c=[];                   % no inequality constraints
    ceq=[outm outn outarma];        % mean, normalization, arma and residual constraints
    gradc = [];
    gradceq = [gradoutm gradoutn gradoutarma] ; % all gradients in one vector
    
else % do not compute gradients if not requested
    for i=1:n
        pexponenti = p.*exp(b(i)+xi(i)*Y') ;
        ypexponenti = Y'.*p.*exp(b(i)+xi(i)*Y') ;
        mmu(i)= sum(ypexponenti) ;
        nnorm(i)= sum(pexponenti) ;
    end
    outm =transpose(mu)-mmu ; % outm is a horizontal vector
    outn = 1 - nnorm ;      % outn is a horiztonal vector
    outarma = transpose(Z)- transpose(Zc);            % arma constraints
    c=[];                   % no inequality constraints
    ceq=[outm outn outarma];        % mean and normalization constraints
end
end