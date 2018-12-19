function [] = plotspglarma(Y, X, link,tests,fitted,phat,sdhat,residualtype)


% four plots are currently available: 
% 1) a time series plot with observed values of the dependent variable, 
%    fixed effects fit (SP-GLM), and GLARMA (SP-GLARMA)fit; 
% 2) an ACF plot of residuals; 
% 3) a plot of residuals against time; 
% 4) the PIT histogram; 
% 
% By default, all 4 plots are provided.
% 
% fitted - fitted values of the SP-GLARMA fit
% phat - phat output from SP-GLARMA 
% sdhat - sdhat output form SP-GLARMA
% residualtype - type of the residual - Pearson / Score (for getting the
% right title) 

switch residualtype 
    case 'Pearson' 
        [~, ~, fitted_fix, ~, ~, ~] = spglm4(Y,X, link, tests); % fitted effects fit
        n = length(Y);
        % ts plot
        subplot(2,2,1)
        h = plot(1:n,Y, 1:n,fitted_fix,1:n,fitted);
        title('Observed vs Fixed vs GLARMA','FontSize',10);
        set(h(1),'LineWidth',1.5,'color','k','LineStyle', '-');
        set(h(2),'LineWidth',1.5,'color','red');
        set(h(2),'LineWidth',1.5,'color','blue');
        ymin = min([Y;fitted_fix;fitted]);
        ymax = max([Y;fitted_fix;fitted]);
        axis([1 n ymin ymax]);
        ylabel('Y','FontSize',10)
        xlabel('Time','FontSize',10)
        legend('Observed','Fixed','GLARMA')
        % acf plot
        subplot(2,2,2)
        resid = (Y-fitted)./sdhat;
        autocorr(resid)
        title('ACF of Pearson Residuals','FontSize',10);
        ylabel('Sample ACF','FontSize',10);
        xlabel('Lag','FontSize',10);
        % resid plot
        subplot(2,2,3)
        plot(1:n,resid)
        axis([1 n min(resid) max(resid)]);
        title('Pearson Residuals','FontSize',10);
        ylabel('Residuals','FontSize',10)
        xlabel('Time','FontSize',10)
        % PIT plot
        subplot(2,2,4)
        SPPIT_hist(Y,phat);
        title('PIT for SP GLARMA','FontSize',10);
        xlabel('Probability Integral Transform','FontSize',10)
        ylabel('Relative Frequency','FontSize',10)
        

    case 'Score'
        [~, ~, fitted_fix, ~, ~, ~] = spglm4(Y,X,link, tests); % fitted effects fit
        n = length(Y);
        % ts plot
        subplot(2,2,1)
        h = plot(1:n,Y, 1:n,fitted_fix,1:n,fitted);
        title('Observed vs Fixed vs GLARMA','FontSize',10);
        set(h(1),'LineWidth',1.5,'color','k','LineStyle', '-');
        set(h(2),'LineWidth',1.5,'color','red');
        set(h(2),'LineWidth',1.5,'color','blue');
        ymin = min([Y;fitted_fix;fitted]);
        ymax = max([Y;fitted_fix;fitted]);
        axis([1 n ymin ymax]);
        ylabel('Y','FontSize',10)
        xlabel('Time','FontSize',10)
        % acf plot
        subplot(2,2,2)
        resid = (Y-fitted)./sdhat;
        autocorr(resid)
        title('ACF of Score Residuals','FontSize',10);
        ylabel('Sample ACF','FontSize',10);
        xlabel('Lag','FontSize',10);
        % resid plot
        subplot(2,2,3)
        plot(1:n,resid)
        axis([1 n min(resid) max(resid)]);
        title('Score Residuals','FontSize',10);
        ylabel('Residuals','FontSize',10)
        xlabel('Time','FontSize',10)
        % PIT plot
        subplot(2,2,4)
        SPPIT_hist(Y,phat);
        title('PIT for SP GLARMA','FontSize',10);
        xlabel('Probability Integral Transform','FontSize',10)
        ylabel('Relative Frequency','FontSize',10)
        
    otherwise
        warning('Unexpected residual type')
end








