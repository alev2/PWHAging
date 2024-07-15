function [Mu] = assemble_Mortality_Matrix(lifeTable,ageMesh)

%This takes the life table as input and returns a diagonal matrix which properly
%accounts for age-dependent mortality

    interp_LifeTable=interp1(lifeTable(:,1),lifeTable(:,2),ageMesh,'linear','extrap');
    interp_LifeTable(isnan(interp_LifeTable))=max(interp_LifeTable);
    
    Mu=diag(-log(1-interp_LifeTable));

end

