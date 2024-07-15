PWHBaseData=readtable('BasicATLASData.csv');

PWHBaseData.Cases=strrep(PWHBaseData.Cases,',','');
PWHBaseData.Cases=cellfun(@str2num,PWHBaseData.Cases);

PWHBaseData.Year=extractBetween(PWHBaseData.Year,1,4);
PWHBaseData.Year=cellfun(@str2num,PWHBaseData.Year);


PWHBaseData.Indicator=strrep(PWHBaseData.Indicator,'HIV ','');
PWHBaseData.Indicator= char(PWHBaseData.Indicator);

for i=1:size(PWHBaseData.Indicator,1)
    PWHBaseData.Indicator(i,:)=replace(PWHBaseData.Indicator(i,:),PWHBaseData.Indicator(i,1),upper(PWHBaseData.Indicator(i,1)));
end



%PWHNewDiagData.Indicator=strrep(PWHNewDiagData.Indicator,'"','');
%PWHDeathData.Year(isnan(PWHDeathData.Year))=2
PWHBaseData=sortrows(PWHBaseData,'Year');

DiagData=PWHBaseData(contains(string(PWHBaseData.Indicator),'Diagnoses')==1,:);
PrevData=PWHBaseData(contains(string(PWHBaseData.Indicator),'Prevalence')==1,:);
DeathData=PWHBaseData(contains(string(PWHBaseData.Indicator),'Deaths')==1,:);

%load('PWHNewDiagData.mat');
mode1=PrevData.Cases(1:end-1)+(DiagData.Cases(1:end-1)-DeathData.Cases(1:end-1));
mode2=PrevData.Cases(1:end-1)+(DiagData.Cases(2:end)-DeathData.Cases(2:end));
PrevUpdate=PrevData.Cases(2:end);

YrUpdate=PrevData.Year(2:end);

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



tp=tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile

plot(YrUpdate,PrevUpdate,'--','LineWidth',lw,'Color',colorMat(1,:));
hold on
plot(YrUpdate,mode1,'--','LineWidth',lw,'Color',colorMat(2,:));
plot(YrUpdate,mode2,'--','LineWidth',lw,'Color',colorMat(3,:));


legend({'Surv data','Mode1','Mode2'},'FontSize',legendfs,'Interpreter','latex','location','northwest','NumColumns',1);
xlabel('Year','FontSize',fs,'Interpreter','latex');
ylabel('Annual HIV Prevalence [Persons]','FontSize',fs,'Interpreter','latex');
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


