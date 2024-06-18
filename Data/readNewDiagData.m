PWHNewDiagData=readtable('PWHNewDiagData.csv');
PWHNewDiagData.Cases=strrep(PWHNewDiagData.Cases,',','');
PWHNewDiagData.Cases=cellfun(@str2num,PWHNewDiagData.Cases);

PWHNewDiagData.AgeGroup=extractBetween(PWHNewDiagData.AgeGroup,1,2);
PWHNewDiagData.AgeGroup=cellfun(@str2num,PWHNewDiagData.AgeGroup);

PWHNewDiagData.Indicator=string(     char(PWHNewDiagData.Indicator));
%PWHNewDiagData.Indicator=strrep(PWHNewDiagData.Indicator,'"','');
PWHNewDiagData.Year(isnan(PWHNewDiagData.Year))=2020;
PWHNewDiagData=sortrows(PWHNewDiagData,'AgeGroup');
PWHNewDiagData=sortrows(PWHNewDiagData,'Year');
%load('PWHNewDiagData.mat');
PWHNewDiagData.SaxphonePhil=zeros(size(PWHNewDiagData,1),1);
PWHNewDiagData.TromboneTony=zeros(size(PWHNewDiagData,1),1);