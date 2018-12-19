function [c, ceq, gradc, gradceq] = constraints4(para, x, y, link)
% computes the mean and normalization constraints
% y is a vertical n-vector, x is a n by p matrix

n = length(y) ;
q = size(x,2) ;
beta = para(1:q) ;              %mean parameters
p = exp(para(q+1:q+n)) ;        %probability vector
b = para(q+n+1: q+2*n) ;        %normalizing constants
theta =para(q+2*n+1:q+3*n) ;    %tilt values

% mean values under beta and test constraints
if  strcmp(link,'id')
    mu = x*beta; 
elseif  strcmp(link,'log')
    mu = exp(x*beta);
elseif  strcmp(link,'inv')
    mu = 1./(x*beta);
elseif  strcmp(link,'logit')
    mu = exp(x*beta)./(1+exp(x*beta)) ; 
end;

% to store mean and norm constraints
% can be written more efficiently without for loops
mmu = zeros(1,n);
nnorm = zeros(1,n) ;

% compute mean, norm constraints and their gradients (if requested)
if nargout > 2
    gradoutm = zeros(q+3*n, n) ; % column i = mean constraint i
    gradoutn = zeros(q+3*n, n) ; % column i = norm constraint i
    
    for i=1:n
        % precompute each of these quantities
        pexponenti = p.*exp(b(i)+theta(i)*y) ;
        ypexponenti = y.*p.*exp(b(i)+theta(i)*y) ;
        y2pexponenti = y.^2.*p.*exp(b(i)+theta(i)*y) ;
        
        % mean and norm constraints
        mmu(i)= sum(ypexponenti) ;
        nnorm(i)= sum(pexponenti) ;
        
        % compute gradients of mean constraint
        if  strcmp(link,'id')
            gradoutm(:,i) = [x(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'log')
            gradoutm(:,i) = [mu(i)*x(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'inv')
            gradoutm(:,i) = [-mu(i)^2*x(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        elseif  strcmp(link,'logit')
            gradoutm(:,i) = [mu(i)*(1-mu(i))*x(i,:)' ; -ypexponenti ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]' ; [zeros(1,i-1) -sum(y2pexponenti) zeros(1, n-i)]'] ;
        end;
        
        % compute gradient of norm constraint
        gradoutn(:,i) = [zeros(q,1) ; -pexponenti ; [zeros(1,i-1) -nnorm(i) zeros(1, n-i)]' ; [zeros(1,i-1) -mmu(i) zeros(1, n-i)]'] ;
    end;
    
    outm = transpose(mu)-mmu ; % outm is a horizontal vector
    outn = 1 - nnorm ;      % outn is a horiztonal vector
    c=[];                   % no inequality constraints
    ceq=[outm outn];        % mean and normalization constraints
    gradc = [];
    gradceq = [gradoutm gradoutn] ; % all gradients in one vector
    
else % do not compute gradients if not requested
    for i=1:n
        pexponenti = p.*exp(b(i)+theta(i)*y) ;
        ypexponenti = y.*p.*exp(b(i)+theta(i)*y) ;
        mmu(i)= sum(ypexponenti) ;
        nnorm(i)= sum(pexponenti) ;
    end;
    outm =transpose(mu)-mmu ; % outm is a horizontal vector
    outn = 1 - nnorm ;      % outn is a horiztonal vector
    c=[];                   % no inequality constraints
    ceq=[outm outn];        % mean and normalization constraints
end;
end