function [mortRatesCopy] = generateInvKF_BarPlot_SubPts(deathsByYearReal,deathsByYearKF,years,plotYears,ageMesh,mortRates,ageMesh2,interpType)
%GENERATEINVKF_BARPLOT Summary of this function goes here
%   Detailed explanation goes here
    lw=3;
    fs=18;
    legendfs=14;
    ms=10;

%   for 4 by 4    
%    posMat=[-1000 0 1500 650];
    
%   for surv deaths    .
    posMat=[0 0 1500 350];

%   for curves
%    posMat=[0 0 1500 500]; legendfs=18;

%for the mortality surface and 
%    posMat=[0 0 800 450];
    
    mortRatesCopy=[];
    
    
    
    for kk=2:size(mortRates,2)
       mortRateCur=interpn(ageMesh2,mortRates(:,kk),ageMesh,interpType);
       mortRatesCopy=[mortRatesCopy mortRateCur];
    end
    
    mortRates=mortRatesCopy;   
    
    
    [~,plotInds]=intersect(years,plotYears);
    
    
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

    tp=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    
%     nexttile
%     
%     bar(years,deathsByYearReal');
%     legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
%     title('Observed mortality','FontSize',fs,'Interpreter','latex');
%     xlabel('Year','FontSize',fs,'Interpreter','latex');
%     ylabel('PWH Deaths [Deaths/year]','FontSize',fs,'Interpreter','latex');
%     ylim([0 9000]);
%     set(gca,'FontSize',fs);
%     ax=gca;
%     ax.YAxis.Exponent=0;
%     
%     nexttile
%     
%     bar(years,deathsByYearKF');
%     title('Simulated mortality','FontSize',fs,'Interpreter','latex');
%     legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
%     xlabel('Year','FontSize',fs,'Interpreter','latex');
%     ylabel('PWH Deaths [Deaths/year]','FontSize',fs,'Interpreter','latex');
%     set(gca,'FontSize',fs);
%     ax=gca;
%     ax.YAxis.Exponent=0;
%     ylim([0 9000]);

    
    
%     nexttile
%     
%     
%    % bottomMedicine=...
%    %     [0.0 0.0 .25;...
%    %      0.0 0.5 1.0;...
%    %      1.0 1.0 1.0;...
%    %      1.0 0.0 0.0;...
%    %      0.5 0.0 0.0];
% 
%      
% %     bottomMedicine=...
% %         [0.0 0.0 .5;...
% %          0.0 0.5 1.0;...
% %          0.0 1.0 0.0;...
% %          %0.5 1.0 0.5;...
% %          %1.0 1.0 1.0;...
% %          0.34 0.66 0.0;...
% %          0.5  0.5 0.0;...
% %          1.0 1.0 0.0;...        
% %          1.0 0.5 0.0;...
% %          1.0 0.0 0.0;...
% %          0.5 0.0 0.0];    
%      
% %    bottomMedicine=...
% %         [0.3010 0.7450 0.9330;...
% %          ...0.0 0.0 1.0;...
% %          ...0.0 0.5 1.0;...
% %          0.0 1.0 0.0;...
% %          %0.5 1.0 0.5;...
% %          %1.0 1.0 1.0;...
% %          %0.34 0.66 0.0;...
% %          ...0.5  0.5 0.0;...
% %          1.0 1.0 0.0;...        
% %          ...%1.0 0.5 0.0;...
% %          1.0 0.0 0.0;...
% %          ...%0.5 0.0 0.0
% %          ];          
%      
%    bottomMedicine=...
%         [
%          0.0 0.0 1.0;...
%          0.0 0.5 1.0;...
%          0.0 1.0 1.0;...
%          0.0 1.0 0.5;...
%          0.0 1.0 0.0;...
%          0.5 1.0 0.0;...
%          1.0 1.0 0.0;...        
%          1.0 0.5 0.0;...        
%          1.0 0.0 0.0;...
%          ];              
%      
%     bottomMedicine2=[];
% 
% 	xs=linspace(0,1,10);
%     for i=1:(size(bottomMedicine,1)-1)
% 
%         for j=1:length(xs)
%             bottomMedicine2=...
%                 [bottomMedicine2;...
%                  (1-xs(j))*bottomMedicine(i,:)+xs(j)*bottomMedicine(i+1,:)];
%         end
%         
%     end
%    
%     
%     surf(ageMesh(1:4:(end-10)),2009:1:2022,1-exp(-mortRates(1:4:(end-10),1:end))')
%     h=gca;
%     xlim([40 85])
%     xlabel('Age','Interpreter','latex','FontSize',fs);    
%     ylim([2009 2022])
%     ylabel('Year','Interpreter','latex','FontSize',fs);
%     zlim([.008 .15]); 
%     zlabel('Annual prob. mortality [-]','Interpreter','latex');
%     set(h,'ZScale','log');
%     c=colorbar('location','southoutside');
%     c.Label.String='Annual prob. mortality [-]';
%     c.Label.Interpreter='latex';
%     c.Label.FontSize=fs-3;
%     colormap(bottomMedicine2);
%     set(gca,'FontSize',fs-3);
%     caxis([.01 .05]);
    

    
    plotStyles={...
     '-',...
     '-.',...
     '--',...
     '-o',...
     ':',...
     '--^',...
     '-.s',...
     ':d'...
     };
    
    lineWidths=[lw-1;lw;lw;lw;lw+3; lw; lw; lw];
    skipFactors=[1;1;1;8;1;8;8;8];
    nexttile
    for kk=1:length(plotInds)

        semilogy(ageMesh(1:skipFactors(kk):end-1),1-exp(-mortRates(1:skipFactors(kk):end-1,plotInds(kk))),plotStyles{kk},'Color',colorMat(kk,:),'LineWidth',lineWidths(kk));
        if(kk==1)
            hold on
        end
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,4)),'--','Color',colorMat(2,:),'LineWidth',lw);
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,6)),'--','Color',colorMat(3,:),'LineWidth',lw);
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,8)),'--','Color',colorMat(4,:),'LineWidth',lw);
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,10)),'--','Color',colorMat(5,:),'LineWidth',lw);    
    end
    
    legend(num2str(years(plotInds)),'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
%    legend({'2010','2012','2014','2016','2018'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
    xlabel('Age','FontSize',fs,'Interpreter','latex');
    ylabel('Annual prob. mortality [-]','FontSize',fs,'Interpreter','latex');
   
    title('Prob. mortality by age','FontSize',fs,'Interpreter','latex');
%   legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',fs,'Interpreter','latex','location','northwest','NumColumns',2);
    set(gca,'FontSize',fs);
    ax=gca;
    ax.YAxis.Exponent=0;
    ylim([1e-5 1]);
%    ylim([1e-3 .3]);
    yticks([.00001 .0001 .001 .01 .1 ]);
    
    xlim([0 101]);
    
    nexttile
    
    for kk=1:length(plotInds)

        plot(ageMesh(1:skipFactors(kk):end-1),1-exp(-mortRates(1:skipFactors(kk):end-1,plotInds(kk))),plotStyles{kk},'Color',colorMat(kk,:),'LineWidth',lineWidths(kk));
        if(kk==1)
            hold on
        end
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,4)),'--','Color',colorMat(2,:),'LineWidth',lw);
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,6)),'--','Color',colorMat(3,:),'LineWidth',lw);
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,8)),'--','Color',colorMat(4,:),'LineWidth',lw);
%        semilogy(ageMesh(1:end-1),1-exp(-mortRates(1:end-1,10)),'--','Color',colorMat(5,:),'LineWidth',lw);    
    end
    
    legend(num2str(years(plotInds)),'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
%    legend({'2010','2012','2014','2016','2018'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',2);
    xlabel('Age','FontSize',fs,'Interpreter','latex');
    ylabel('Annual prob. mortality [-]','FontSize',fs,'Interpreter','latex');
   
    title('Prob. mortality by age','FontSize',fs,'Interpreter','latex');
%   legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',fs,'Interpreter','latex','location','northwest','NumColumns',2);
    set(gca,'FontSize',fs);
    ax=gca;
    ax.YAxis.Exponent=0;
%    ylim([1e-4 1]);
    ylim([.0085 .055]);
    yticks([.0075 .01 .025 .05 .055 ]);
    
    xlim([40 80]);
    



    set(gcf,'Position',posMat);   

end

