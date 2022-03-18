%inputs:fig is which figures it creates
%       ylimit is 2x2 array, rows describing limits of grpahs' y axes

function createFnMProfiles2n3(avFMRPM,SEM,avFMRPMCorr,SEMCorr, fig, ylimit)
figure(fig(1))
%title('T-motor Coplanar Force Profiles') <------Add title in powerpoint
axistitles=["Fx (N)", "Fy (N)","Fz (N)"];
for i=1:3 %loop through Fx, Fy, Fz
    subplot(3,1,i)
    hold on;
    grid on;
    %errorbar(avFMRPM(:,14),avFMRPM(:,i+6),SEM(:,i+6),'--xr'  ) %plots Fx vs RPM for uncorrected
    errorbar(avFMRPMCorr(:,14),avFMRPMCorr(:,i+6),SEMCorr(:,i+6),'--xk'  ) %plots for corrected
    %legend({'Uncorrected','Corrected'},'Location','northwest')
    xlim([2500,5000])
    ylim(ylimit(1,:))
    xlabel('RPM')
    ylabel(axistitles(i))
    hold off;
end

figure(fig(2))
%title('T-motor Coplanar Moment Profiles') <-----Add title in powerpoint
axistitles=["Mx (Nm)","My (Nm)","Mz (Nm)"];
for i=1:3 %loop through Mx, My, Mz
    subplot(3,1,i)
    hold on;
    grid on;
    %errorbar(avFMRPM(:,14),avFMRPM(:,i+9),SEM(:,i+9),'-xr'  ) %plots Fx vs RPM for uncorrected
    errorbar(avFMRPMCorr(:,14),avFMRPMCorr(:,i+9),SEMCorr(:,i+9),'--xk'  ) %plots for corrected
    %legend({'Uncorrected','Corrected'},'Location','southwest')
    xlim([2500,5000])
    ylim(ylimit(2,:))
    xlabel('RPM')
    ylabel(axistitles(i))
    hold off;
end

end