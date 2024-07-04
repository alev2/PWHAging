PWHNewDiagData_55andUp=readtable('PWHNewDiagData_55andUp.csv');
PWHNewDiagData_55andUp.Cases=strrep(PWHNewDiagData_55andUp.Cases,',','');
PWHNewDiagData_55andUp.Cases=cellfun(@str2num,PWHNewDiagData_55andUp.Cases);

PWHNewDiagData_55andUp.AgeGroup=extractBetween(PWHNewDiagData_55andUp.AgeGroup,1,2);
PWHNewDiagData_55andUp.AgeGroup=cellfun(@str2num,PWHNewDiagData_55andUp.AgeGroup);

PWHNewDiagData_55andUp.Indicator=string(     char(PWHNewDiagData_55andUp.Indicator));
%PWHNewDiagData.Indicator=strrep(PWHNewDiagData.Indicator,'"','');
%PWHNewDiagData_55andUp.Year(isnan(PWHNewDiagData_55andUp.Year))=2020;
PWHNewDiagData_55andUp.Year=extractBetween(PWHNewDiagData_55andUp.Year,1,4);
PWHNewDiagData_55andUp.Year=cellfun(@str2num,PWHNewDiagData_55andUp.Year);
PWHNewDiagData_55andUp.RatePer100000=cellfun(@str2num,PWHNewDiagData_55andUp.RatePer100000);
PWHNewDiagData_55andUp=sortrows(PWHNewDiagData_55andUp,'AgeGroup');
PWHNewDiagData_55andUp=sortrows(PWHNewDiagData_55andUp,'Year');
%load('PWHNewDiagData.mat');