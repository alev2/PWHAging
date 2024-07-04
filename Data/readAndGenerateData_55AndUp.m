

PWHPrevData55Up=readtable('PWHPrevData_55andUp.csv');
PWHPrevData55Up.Cases=strrep(PWHPrevData55Up.Cases,',','');
PWHPrevData55Up.Cases=cellfun(@str2num,PWHPrevData55Up.Cases);

PWHPrevData55Up.AgeGroup=extractBetween(PWHPrevData55Up.AgeGroup,1,2);
PWHPrevData55Up.AgeGroup=cellfun(@str2num,PWHPrevData55Up.AgeGroup);

PWHPrevData55Up.Indicator=string(     char(PWHPrevData55Up.Indicator));
%PWHPrevData.Indicator=strrep(PWHPrevData.Indicator,'"','');
PWHPrevData55Up.Year=extractBetween(PWHPrevData55Up.Year,1,4);
PWHPrevData55Up.Year=cellfun(@str2num,PWHPrevData55Up.Year);
PWHPrevData55Up=sortrows(PWHPrevData55Up,'Year');
PWHPrevData55Up=sortrows(PWHPrevData55Up,'AgeGroup');
%load('PWHPrevData.mat');
PWHPrevData55Up.SaxphonePhil=zeros(size(PWHPrevData55Up,1),1);
PWHPrevData55Up.TromboneTony=zeros(size(PWHPrevData55Up,1),1);

yrStart=2008;
yrArray=(2008:1:2023);

[~,startInd]=find(yrArray==yrStart);

totInc=[];
incPct=[];
barPlotHell=[];

for i=2008:2023
   
    incTmp=sum(PWHPrevData55Up.Cases(PWHPrevData55Up.Year==i));

    PWHPrevData55Up.SaxphonePhil(PWHPrevData55Up.Year==i) =...
        PWHPrevData55Up.Cases(PWHPrevData55Up.Year==i)/incTmp;

    PWHPrevData55Up.TromboneTony(PWHPrevData55Up.Year==i)=incTmp;

    barPlotHell=[barPlotHell; PWHPrevData55Up.SaxphonePhil(PWHPrevData55Up.Year==i)'];
    
end




bar(100*barPlotHell(startInd:end,:));

set(gca,'FontSize',fs);

ax=ancestor(gca,'Axes');
ax.YAxis.Exponent=0;
title('Age distribution PWH, 2012-22',...
        'Interpreter','latex','FontSize',fs);

    
    
xticklabels(arrayfun(@(a)num2str(a),yrArray(startInd:end),'uni',0) );

legend({'13-24','25-34','35-44','45-54','55-64','65+'},'FontSize',fs,'Interpreter','latex','location','northwest','NumColumns',2);

ylabel('\% of PWH [-]', 'Interpreter','latex','FontSize',fs);
xlabel('Year', 'Interpreter','latex','FontSize',fs);


