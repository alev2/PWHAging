function [] = plotTimeSeries(solutionMatrix,timeOutput,ageMesh)
%
% 
% 
lw=3;
fs=15;

 
for i=1:size(solutionMatrix,2) 

    plot(ageMesh,solutionMatrix(:,i),'LineWidth',lw);

    ax=ancestor(gca,'Axes');
    ax.YAxis.Exponent=0;
    set(gca,'FontSize',fs);
    axis([0 100  0 35000]);

    title(strcat('PWH Population distribution: ',num2str(floor(timeOutput(i)))), 'FontSize',fs,'Interpreter','latex');

    xlabel('Age','FontSize',fs,'Interpreter','latex');
    ylabel('Num. PWH','FontSize',fs,'Interpreter','latex');

    input('Press enter to continue');

end

