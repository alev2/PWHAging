function [outputArg1,outputArg2] = convergencePlots(convergenceHistories,inclYrs)
%CONVERGENCEPLOTS Summary of this function goes here
%   Detailed explanation goes here

    if(nargin==1)
        inclYrs=convergenceHistories(1,:)';
    end
    lw=3;
    fs=18;
    legendfs=14;
    ms=10;
    
    posMat=[-1000 0 1500 650];
    

    colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184;
    .1235 .078 .484;
    .8235 .3278 .384;
    .5235 .7278 .184;
    .11 .66 .55;
    .51 .26 .85;
    .51 .96 .05;
    .51 .556 .15;
    .151 .356 .75;
    0 0 0    
    ];

    [~,subsetInds]=intersect(convergenceHistories(1,:)',inclYrs);

    convYears=num2str(convergenceHistories(1,subsetInds)');  
    
    
    
    numIts=(1:1:size(convergenceHistories(2:end,:),1))';    
    convergenceHistories=convergenceHistories(2:end,:);
    
    tp=tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    nexttile
    size(numIts)
    
    
    semilogy(numIts,convergenceHistories(:,subsetInds),':','LineWidth',lw);
    
    title('Convergence by year','FontSize',fs,'Interpreter','latex');
    legend(convYears,'FontSize',legendfs,'Interpreter','latex','location','northeast','NumColumns',3);
    
    xlabel('Iteration','FontSize',fs,'Interpreter','latex');
    ylabel('Conv. Functional','FontSize',fs,'Interpreter','latex');

    set(gca,'FontSize',fs);
    ax=gca;
    ax.YAxis.Exponent=0;

end

