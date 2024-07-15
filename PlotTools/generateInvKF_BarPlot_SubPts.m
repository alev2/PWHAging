function [] = generateInvKF_BarPlot_SubPts(deathsByYearReal,deathsByYearKF,years,ageMesh,mortRates,ageMesh2,interpType)
%GENERATEINVKF_BARPLOT Summary of this function goes here
%   Detailed explanation goes here
    lw=3;
    fs=18;
    legendfs=14;
    ms=10;
    
    posMat=[-1000 0 1500 650];
    
    mortRatesCopy=[];
    for kk=1:size(mortRates,2)
       mortRateCur=interpn(ageMesh2,mortRates(:,kk),ageMesh,interpType);
       mortRatesCopy=[mortRatesCopy mortRateCur];
    end
    mortRates=mortRatesCopy;   
    
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

    tp=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile
    
    bar(years,deathsByYearReal');
    legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
    title('Observed mortality','FontSize',fs,'Interpreter','latex');
    xlabel('Year','FontSize',fs,'Interpreter','latex');
    ylabel('PWH Deaths [Deaths/year]','FontSize',fs,'Interpreter','latex');
    ylim([0 9000]);
    set(gca,'FontSize',fs);
    ax=gca;
    ax.YAxis.Exponent=0;
    
    nexttile
    
    bar(years,deathsByYearKF');
    title('Simulated mortality','FontSize',fs,'Interpreter','latex');
    legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
    xlabel('Year','FontSize',fs,'Interpreter','latex');
    ylabel('PWH Deaths [Deaths/year]','FontSize',fs,'Interpreter','latex');
    set(gca,'FontSize',fs);
    ax=gca;
    ax.YAxis.Exponent=0;
    ylim([0 9000]);

    nexttile
    surf(ageMesh(1:4:(end-10)),2009:1:2022,1-exp(-mortRates(1:4:(end-10),2:end))')
    h=gca;
    xlim([40 70])
    set(h,'ZScale','log');
    colorbar;
    colormap('hsv');
    set(gca,'FontSize',fs-3);
    caxis([.01 .1]);
    set(gcf,'Position',posMat);   
    
    nexttile

    semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,2)),'--','Color',colorMat(1,:),'LineWidth',lw);
    hold on
    semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,4)),'--','Color',colorMat(2,:),'LineWidth',lw);
    semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,6)),'--','Color',colorMat(3,:),'LineWidth',lw);
    semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,8)),'--','Color',colorMat(4,:),'LineWidth',lw);
    semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,10)),'--','Color',colorMat(5,:),'LineWidth',lw);
    %semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,12)),'--','Color',colorMat(6,:),'LineWidth',lw);
    %semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,14)),'--','Color',colorMat(7,:),'LineWidth',lw);
    legend({'2010','2012','2014','2016','2018'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
    xlabel('Age','FontSize',fs,'Interpreter','latex');
    ylabel('Annual prob. mortality [-]','FontSize',fs,'Interpreter','latex');
   
    title('Prob. mortality by age','FontSize',fs,'Interpreter','latex');
%   legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',fs,'Interpreter','latex','location','northwest','NumColumns',2);
    set(gca,'FontSize',fs);
    ax=gca;
    ax.YAxis.Exponent=0;
    ylim([3.15e-4 1]);
    


end

