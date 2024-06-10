


load('PWHPrevData.mat');
PWHPrevData.SaxphonePhil=zeros(size(PWHPrevData,1),1);
PWHPrevData.TromboneTony=zeros(size(PWHPrevData,1),1);

yrStart=2012;
yrArray=(2008:1:2021);

[~,startInd]=find(yrArray==yrStart);

totInc=[];
incPct=[];
barPlotHell=[];

for i=2008:2021
   
    incTmp=sum(PWHPrevData.Cases(PWHPrevData.Year==i));

    PWHPrevData.SaxphonePhil(PWHPrevData.Year==i) =...
        PWHPrevData.Cases(PWHPrevData.Year==i)/incTmp;

    PWHPrevData.TromboneTony(PWHPrevData.Year==i)=incTmp;

    barPlotHell=[barPlotHell; PWHPrevData.SaxphonePhil(PWHPrevData.Year==i)'];
    
end




bar(100*barPlotHell(startInd:end,:));

set(gca,'FontSize',fs);

ax=ancestor(gca,'Axes');
ax.YAxis.Exponent=0;
title('Age distribution PWH, 2012-22',...
        'Interpreter','latex','FontSize',fs);

    
    
xticklabels(arrayfun(@(a)num2str(a),yrArray(startInd:end),'uni',0) );

legend({'13-24','25-34','35-44','45-54','55+'},'FontSize',fs,'Interpreter','latex','location','northwest','NumColumns',2);

ylabel('\% of PWH [-]', 'Interpreter','latex','FontSize',fs);
xlabel('Year', 'Interpreter','latex','FontSize',fs);


