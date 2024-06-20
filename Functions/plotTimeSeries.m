function [totalOutput] = plotTimeSeries(solutionMatrix,timeOutput,ageMesh)
%
% 
% 
lw=3;
fs=15;

totalOutput=[]; 
for i=1:size(solutionMatrix,2) 

    plot(ageMesh,solutionMatrix(:,i),'LineWidth',lw);

    totalOutput=[totalOutput;trapz(ageMesh,solutionMatrix(:,i))];

    ax=ancestor(gca,'Axes');
    ax.YAxis.Exponent=0;
    set(gca,'FontSize',fs);
    axis([ageMesh(1) ageMesh(end)  0 40000]);

    title(strcat('PWH Population distribution: ',num2str(floor(timeOutput(i)))), 'FontSize',fs,'Interpreter','latex');

    xlabel('Age','FontSize',fs,'Interpreter','latex');
    ylabel('Num. PWH','FontSize',fs,'Interpreter','latex');

    input('Press enter to continue');

end

