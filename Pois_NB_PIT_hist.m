function [] = Pois_NB_PIT_hist(PIT,Type)

% PIT = PIT from glarmaPIT(glarmamodnb)$PIT
% Type = 'Pois' or 'NegBin' 
% Incorproate Type as part of the title of the histogram and it does
% nothing else.

bin = 10;
n = size(PIT,2);
dummy = 0:(1/bin):1;
height = diff(PIT(:,n))*bin;
if (max(height) > 2)
        yupper = max(height) + 1/(bin/2);
else
        yupper = 2;
end

bar(dummy(1:bin)+0.05,height);

switch Type
    case 'Pois' 
        title('PIT for GLARMA (Poisson) ','FontSize',30);
        %title('Poisson ','FontSize',30);
    case 'NegBin'
        title('PIT for GLARMA (Negative Binomial) ','FontSize',30);
        %title('Negative Binomial ','FontSize',30);
    otherwise
        warning('Unexpected Type')
end


line('XData', [0 1], 'YData', [1 1], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
xlabel('Probability Integral Transform','FontSize',20)
ylabel('Relative Frequency','FontSize',20)
set(gca,'fontsize',20)
axis([0 1 0 yupper]);

end