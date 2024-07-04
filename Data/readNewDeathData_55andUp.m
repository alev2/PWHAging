PWHDeathData_55andUp=readtable('PWHDeathData_55andUp.csv');
PWHDeathData_55andUp.Cases=strrep(PWHDeathData_55andUp.Cases,',','');
PWHDeathData_55andUp.Cases=cellfun(@str2num,PWHDeathData_55andUp.Cases);

PWHDeathData_55andUp.AgeGroup=extractBetween(PWHDeathData_55andUp.AgeGroup,1,2);
PWHDeathData_55andUp.AgeGroup=cellfun(@str2num,PWHDeathData_55andUp.AgeGroup);

PWHDeathData_55andUp.Indicator=string(     char(PWHDeathData_55andUp.Indicator));
%PWHNewDiagData.Indicator=strrep(PWHNewDiagData.Indicator,'"','');
PWHDeathData_55andUp.Year=extractBetween(PWHDeathData_55andUp.Year,1,4);
PWHDeathData_55andUp.Year=cellfun(@str2num,PWHDeathData_55andUp.Year);
%PWHDeathData.Year(isnan(PWHDeathData.Year))=2020;
PWHDeathData_55andUp=sortrows(PWHDeathData_55andUp,'AgeGroup');
PWHDeathData_55andUp=sortrows(PWHDeathData_55andUp,'Year');
%load('PWHNewDiagData.mat');
