PWHDeathData=readtable('PWHDeathData.csv');
PWHDeathData.Cases=strrep(PWHDeathData.Cases,',','');
PWHDeathData.Cases=cellfun(@str2num,PWHDeathData.Cases);

PWHDeathData.AgeGroup=extractBetween(PWHDeathData.AgeGroup,1,2);
PWHDeathData.AgeGroup=cellfun(@str2num,PWHDeathData.AgeGroup);

PWHDeathData.Indicator=string(     char(PWHDeathData.Indicator));
%PWHNewDiagData.Indicator=strrep(PWHNewDiagData.Indicator,'"','');
PWHDeathData.Year=extractBetween(PWHDeathData.Year,1,4);
PWHDeathData.Year=cellfun(@str2num,PWHDeathData.Year);
%PWHDeathData.Year(isnan(PWHDeathData.Year))=2020;
PWHDeathData=sortrows(PWHDeathData,'AgeGroup');
PWHDeathData=sortrows(PWHDeathData,'Year');
%load('PWHNewDiagData.mat');
