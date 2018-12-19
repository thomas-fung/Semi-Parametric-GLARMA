%% 
% Examples to test the spglarma function
% associated with Huang & Fung (2017)

% associated files: 
% spglarmapearson, spglarmascore, loglikglarma, constraintsglarmapearson,
% constraintsglarmascore, plotspglarma, SPPIT_hist, spglm4, loglik4,
% constraints4. 

%% Example: Polio data of Davis, Dumsmuir and Wang (2000)
data = csvread('./data/polio.csv',1,2);
Y = data(:,1);
X = data(:,2:7);
[delta, maxloglik, fitted, iter, phat, sdhat] = spglarma(Y,X,[],[1,2,5],'log',[],'Pearson');
%%%%%%%%%%%%%%%%%%%%%%%%
% input (in order): 
% Y: response vector
% X: design matrix. Should have a column of 1
% phiLags: AR type lags; no AR lags specified so entered as [] 
% thetaLags: MA type lags; Lags 1,2 & 5 specified 
% 'log': type of link being used; only 'log' link is fully functional 
% 'test'; the test argument allows you to test hypothesis; entered as [] 
% at the moment as no test is being used
% 'Pearson': 'Pearson' or 'Score' type residuals to be used 
%%%%%%%%%%%%%%%%%%%%%%%%
% output (in order):
% delta: estimated coefficients
% maxloglik: maximized log-likelihood
% fitted: fitted mu for each response
% iter: no of iterations till convergence
% phat: tilted density for each response 
% sdhat: fitted sd for each response

% To produce the stair-case graph for Polio data 
% Poisson and Neg Bin estimates are obtained from the glarma package in R
mean1 = fitted(36) ;
phat1 = phat(36,:) ;
[sortedy , ordery] = sort(Y) ;
Fhat1 = cumsum(phat1(ordery));
% 
stairs([sortedy;15;16], [Fhat1,1,1]) ;
hold on; 
stairs(0:16, poisscdf(0:16,   5.939786), 'col','r') ;
mean1 =  6.730633; 
stairs(0:16, nbincdf(0:16, 1.588,1.588./(1.588+mean1)), 'col','g') ;
hold off ;

% To calculate the eq.se for the Polio data set
data = csvread('./data/polio.csv',1,2);
Y = data(:,1);
X = data(:,2:7);
tic
[delta1, maxloglik1] = spglarma(Y,X,[],[1,2,5],'log',[],'Pearson');
[~, maxloglik2] = spglarma(Y,X,[],[1,2,5],'log',...
    [0, NaN*ones(1,8)],'Pearson');
[~, maxloglik3] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,1),0, NaN*ones(1,7)],'Pearson');
[~, maxloglik4] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,2),0, NaN*ones(1,6)],'Pearson');
[~, maxloglik5] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,3),0, NaN*ones(1,5)],'Pearson');
[~, maxloglik6] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,4),0, NaN*ones(1,4)],'Pearson');
[~, maxloglik7] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,5),0, NaN*ones(1,3)],'Pearson');
[~, maxloglik8] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,6),0, NaN*ones(1,2)],'Pearson');
[~, maxloglik9] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,7),0, NaN*ones(1,1)],'Pearson');
[~, maxloglik10] = spglarma(Y,X,[],[1,2,5],'log',...
    [NaN*ones(1,8),0],'Pearson');

equivsex1 = abs(delta1(1))./sqrt(2*(maxloglik1 - maxloglik2))
equivsex2 = abs(delta1(2))./sqrt(2*(maxloglik1 - maxloglik3))
equivsex3 = abs(delta1(3))./sqrt(2*(maxloglik1 - maxloglik4))
equivsex4 = abs(delta1(4))./sqrt(2*(maxloglik1 - maxloglik5))
equivsex5 = abs(delta1(5))./sqrt(2*(maxloglik1 - maxloglik6))
equivsex6 = abs(delta1(6))./sqrt(2*(maxloglik1 - maxloglik7))
equivsex7 = abs(delta1(7))./sqrt(2*(maxloglik1 - maxloglik8))
equivsex8 = abs(delta1(8))./sqrt(2*(maxloglik1 - maxloglik9))
equivsex9 = abs(delta1(9))./sqrt(2*(maxloglik1 - maxloglik10))
round([equivsex1, equivsex2, equivsex3, equivsex4, equivsex5, equivsex6, equivsex7, ...
    equivsex8, equivsex9],3)
toc

% To produce the PIT graph from the Polio data
data = csvread('./data/polio.csv',1,2);
Y = data(:,1);
X = data(:,2:7);
tic
[~, ~, ~, ~, phat] = spglarma(Y,X,[],[1,2,5],'log',[],'Pearson');
toc
PIT1 = dlmread('./results/PIT_Pois_Polio.txt',' ');
PIT2 = dlmread('./results/PIT_NB_Polio.txt',' ');
SPPIT_hist(Y,phat);
title('Semiparametric','FontSize',30)
xlabel('')
Pois_NB_PIT_hist(PIT1,'Pois');
title('Poisson','FontSize',30)
xlabel('')
Pois_NB_PIT_hist(PIT2,'NegBin');
title('Negative binomial','FontSize',30)
xlabel('')

% Calculate the Scores for Polio data

% The function is saved as scoringsp.m 
% usage [logarithmic quadratic spherical rankprob dawseb normsq sqerror] 
% = scoringsp(y,phat, fitted, sdhat)
data = csvread('./data/polio.csv',1,2);
Y = data(:,1);
X = data(:,2:7);
[~, ~, fitted, ~, phat,sdhat] = spglarma(Y,X,[],[1,2,5],'log',[],'Pearson');
plotspglarma(Y,X,'log',[],fitted,phat,sdhat,'Pearson');
results = scoringsp(Y,phat,fitted,sdhat);
results

% To compare the stair case & PIT graphs between SP, Poi, Neg Bin for the Polio data 
data = csvread('./data/polio.csv',1,2);
Y = data(:,1);
X = data(:,2:7);
[delta, maxloglik, fitted1, iter, phat, sdhat] = spglarma(Y,X,[],[1,2,5],'log',[],'Pearson');
fittedpoi = dlmread('./results/fitted_poi_for_polio.txt',' ');
fittednb = dlmread('./results/fitted_nb_for_polio.txt',' ');
fittednbalpha = 1.588; 

figure(1)
subplot(2,3,4)
ind = 12;
phat1 = phat(ind,:) ;
[sortedy , ordery] = sort(Y) ;
Fhat1 = cumsum(phat1(ordery));
stairs([sortedy;15;16], [Fhat1,1,1],'col','k','LineWidth',2) ;
hold on; 
stairs(0:16, poisscdf(0:16, fittedpoi(ind)), '--','col','k',...
    'LineWidth',2) ;
stairs(0:16, nbincdf(0:16, fittednbalpha,...
    fittednbalpha/(fittednbalpha+fittednb(ind))), ':','col','k',...
    'LineWidth',2) ;
hold off ;
set(gca,'fontsize',20)
axis([0 16 0 1.1])
title('(d) Observation 12','FontSize',20)
xlabel('y','FontSize',20)
ylabel('P(Y\leq y)','FontSize',20)
legend('SP','Poisson','Neg Bin','location','southeast');

subplot(2,3,5)

ind = 36;
phat1 = phat(ind,:) ;
[sortedy , ordery] = sort(Y) ;
Fhat1 = cumsum(phat1(ordery));
stairs([sortedy;15;16], [Fhat1,1,1],'col','k','LineWidth',2) ;
hold on; 
stairs(0:16, poisscdf(0:16, fittedpoi(ind)), '--','col','k',...
    'LineWidth',2) ;
stairs(0:16, nbincdf(0:16, fittednbalpha,...
    fittednbalpha/(fittednbalpha+fittednb(ind))), ':','col','k',...
    'LineWidth',2) ;
hold off ;
set(gca,'fontsize',20)
axis([0 16 0 1.1])
title('(e) Observation 36','FontSize',20)
xlabel('y','FontSize',20)
ylabel('P(Y\leq y)','FontSize',20)
legend('SP','Poisson','Neg Bin','location','southeast');
subplot(2,3,6)
ind = 121;
phat1 = phat(ind,:) ;
[sortedy , ordery] = sort(Y) ;
Fhat1 = cumsum(phat1(ordery));
stairs([sortedy;15;16], [Fhat1,1,1],'col','k','LineWidth',2) ;
hold on; 
stairs(0:16, poisscdf(0:16, fittedpoi(ind)), '--','col','k',...
    'LineWidth',2) ;
stairs(0:16, nbincdf(0:16, fittednbalpha,...
    fittednbalpha/(fittednbalpha+fittednb(ind))), ':','col','k',...
    'LineWidth',2) ;
hold off ;
set(gca,'fontsize',20)
axis([0 16 0 1.1])
title('(f) Observation 121','FontSize',20)
xlabel('y','FontSize',20)
ylabel('P(Y\leq y)','FontSize',20)
legend('SP','Poisson','Neg Bin','location','southeast');



PIT1 = dlmread('./results/PIT_Pois_Polio.txt',' ');
PIT2 = dlmread('./results/PIT_NB_Polio.txt',' ');

Pois_NB_PIT_hist(PIT1,'Pois');
title('Poisson','FontSize',20)
xlabel('')
axis([0 1 0 1.57])

Pois_NB_PIT_hist(PIT2,'NegBin');
title('Negative binomial','FontSize',20)
xlabel('')
axis([0 1 0 1.57])

SPPIT_hist(Y,phat);
title('Semiparametric','FontSize',20)
xlabel('')
axis([0 1 0 1.57])
