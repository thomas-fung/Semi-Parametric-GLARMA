function [c, ceq, gradc, gradceq] = constraintsglarmapearson(para, Y, X, phi_lags, theta_lags, link)
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
    mu = 1./(X*beta + Z);
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
        ypexponenti = Y'.*p.*exp(b(i)+xi(i)*Y') ;
        mmu(i)= sum(ypexponenti);
        Varhat(i) = phat(i,:)*((Y-mmu(i)).^2) ;
    end
sdhat = sqrt(Varhat) ;

e = (Y-mmu')./sdhat;

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
    
    gradoutm = zeros(rsq+4*n, n) ; % column i = mean constraint i
    gradoutZc = zeros(rsq+4*n, n); % column i = recursive arma constraint i
    gradoutZctime = [zeros(rsq+4*n, msq) gradoutZc]; %arma constraint with lags
    gradoute = zeros(rsq+4*n,n); % column i = residual i
    gradoutetime = [zeros(rsq+4*n, msq) gradoute]; %residuals with lags
    gradoutn = zeros(rsq+4*n, n) ; % column i = norm constraint i
    gradouta = zeros(rsq+4*n, n) ; % column i = arma constraint i
   
    if s>0
        temptime = zeros(phi_lags(s),1); 
        temptime(phi_lags(s)+1-phi_lags(1:s)) = phi(1:s);
        % to put the coef in reverse order and with the right amounts of 
        % gap
        temptime = [zeros(n-phi_lags(s)-1,1); temptime; zeros(n,1)];
        % add enough zeros before and after, no matter which chuck we are
        % going to select later
        [sizetemptime, ~]= size(temptime);
    end

    for i = 1:n
        
        % precompute each of these quantities
        pexponenti = p.*exp(b(i)+xi(i)*Y') ;
        
        ypexponenti = Y'.*p.*exp(b(i)+xi(i)*Y') ;
        y2pexponenti = Y'.^2.*p.*exp(b(i)+xi(i)*Y');
        y3pexponenti = Y'.^3.*p.*exp(b(i)+xi(i)*Y');
        mmu(i)= sum(ypexponenti);
        sumy2pexponenti = sum(y2pexponenti);
        sumy3pexponenti = sum(y3pexponenti);
        
        gradoute(1:rsq+n,i) = zeros(rsq+n,1);
        gradoute(rsq+n+1:rsq+2*n,i) = -ypexponenti'./sdhat(i)-(0.5.*(Y(i)-mmu(i)).*(y2pexponenti'-2.*mmu(i).*ypexponenti'))./(sdhat(i).^3);
        gradoute(rsq+2*n+1:rsq+3*n,i) = [zeros(i-1,1); -mmu(i)./sdhat(i)- 0.5.*(Y(i)-mmu(i))./(sdhat(i).^3).*(sumy2pexponenti-2.*mmu(i).^2);zeros(n-i,1)];
        gradoute(rsq+3*n+1:rsq+4*n,i) = [zeros(i-1,1); -sumy2pexponenti./sdhat(i)-0.5.*(Y(i)-mmu(i))./(sdhat(i).^3).*(sumy3pexponenti-2*mmu(i).*sumy2pexponenti);zeros(n-i,1)]; 
        gradoutetime(:,i+msq) =  gradoute(:,i);
    
    
        % Zc w.r.t. beta 
        gradoutZc(1:r,i) = zeros(r,1);
        
        % Zc w.r.t. AR and MA components
        if s>0 && q>0  % both AR, MAparameters and 
            for j = 1:s
                gradoutZc(r+j,i) = Ztime(i+msq-phi_lags(j))+etime(i+msq-phi_lags(j));
            end
            for j = 1:q
                gradoutZc(r+s+j,i) = etime(i+msq-theta_lags(j));
            end
        elseif s>0 && q<=0 % AR parameters/ no MA parameters
            for j = 1:s
                gradoutZc(r+j,i) = Ztime(i+msq-phi_lags(j))+etime(i+msq-phi_lags(j));
            end
        else % MA parameters/ no AR parameters
            for j = 1:q
                gradoutZc(r+j,i) = etime(i+msq-theta_lags(j));
            end 
        end
        
        % Zc_t w.r.t. Z_i 
        % they are 0 if there is no AR component
        % if there are some: 
        if s>0
            gradoutZc((rsq+1):(rsq+n),i) = temptime(sizetemptime-n+2-i:sizetemptime+1-i);
        end
        % Z_t w.r.t. p_j, b_i and psi_j (though e_t)
        
        if s>0
            for k = 1:s
                gradoutZc(rsq+n+1:rsq+4*n,i) = gradoutZc(rsq+n+1:rsq+4*n,i) + phi(k).*gradoutetime(rsq+n+1:rsq+4*n,i+msq-phi_lags(k));
            end
        end
        if q>0
            for k = 1:q
                gradoutZc(rsq+n+1:rsq+4*n,i) = gradoutZc(rsq+n+1:rsq+4*n,i) + theta(k).*gradoutetime(rsq+n+1:rsq+4*n,i+msq-theta_lags(k));
            end
        end
        
        
        gradoutZctime(:,i+msq) = gradoutZc(:,i);
    
         
        % mean and norm constraints
        mmu(i)= sum(ypexponenti) ;
        nnorm(i)= sum(pexponenti) ;
        
        % compute gradients of mean constraint
        if  strcmp(link,'id')
            gradoutm(:,i) = [X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'log')
            gradoutm(1:r,i) = mu(i)*X(i,:)';
            gradoutm(r+1:rsq,i) = zeros(s+q,1);
%             if s>0 && q>0
%                 for j = 1:s
%                     gradoutm(r+j,i) = mu(i).*(Ztime(i+msq-phi_lags(j))+ etime(i+msq-phi_lags(j)));
%                 end
%                 for j = 1:q
%                     gradoutm(r+s+j,i) = mu(i).*etime(i+msq-theta_lags(j));
%                 end
%             elseif s>0 && q<0
%                 for j = 1:s
%                     gradoutm(r+j,i) = mu(i).*(Ztime(i+msq-phi_lags(j))+ etime(i+msq-phi_lags(j)));
%                 end
%             else
%                 for j = 1:q
%                     gradoutm(r+j,i) = mu(i).*etime(i+msq-theta_lags(j));
%                 end
%             end
            gradoutm(rsq+1:rsq+n,i) = [zeros(i-1,1); mu(i); zeros(n-i,1)];
            gradoutm(rsq+n+1:rsq+4*n,i) = [-ypexponenti' ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sumy2pexponenti zeros(1, n-i)]'] ;
            
        elseif  strcmp(link,'inv')
            gradoutm(:,i) = [-mu(i)^2*X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'logit')
            gradoutm(:,i) = [mu(i)*(1-mu(i))*X(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        end
        
        % compute gradient of norm constraint
        gradoutn(:,i) = [zeros(rsq+n,1) ; -pexponenti' ; [zeros(1,i-1) -nnorm(i) zeros(1, n-i)]' ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]'] ;

	% compute gradient of the arma constraint
        gradouta(:,i) = [zeros(rsq,1); zeros(i-1,1); 1; zeros(n-i,1);zeros(3*n,1)] - gradoutZc(:,i);
    end
    
    outm = transpose(mu)-mmu ; % outm is a horizontal vector
    outn = 1 - nnorm ;      % outn is a horiztonal vector
    outarma = transpose(Z)- transpose(Zc);            % arma constraints
    c=[];                   % no inequality constraints
    ceq=[outm outn outarma];        % mean, normalization, arma and residual constraints
    gradc = [];
    gradceq = [gradoutm gradoutn gradouta] ; % all gradients in one vector
    
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