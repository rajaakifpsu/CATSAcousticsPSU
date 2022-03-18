%Purpose: to create profile(s) 1 for Tmotor case (See ______ document).
%inputs:self explanatory except for fig. Fig is the figure numbers to use

function createFnMProfiles1(avFMRPM,SEM,avFMRPMCorr,SEMCorr, fig)

figure(fig(1))
grid on;
hold on;
errorbar(avFMRPM(:,14),avFMRPM(:,9),SEM(:,9),'--xr'  ) %plots Fz vs RPM for uncorrected
errorbar(avFMRPMCorr(:,14),avFMRPMCorr(:,9),SEMCorr(:,9),'--xk'  ) %plots for corrected
legend('Uncorrected','Corrected','Location','northwest')
xlim([1500,5000])
title('T-motor Coplanar Thrust Profile')
xlabel('RPM')
ylabel('Thrust (N)')
hold off;

figure (fig(2))
grid on;
hold on;
errorbar(avFMRPM(:,14),avFMRPM(:,12),SEM(:,12),'--xr'  ) %plots Mz vs RPM for uncorrected
errorbar(avFMRPMCorr(:,14),avFMRPMCorr(:,12),SEMCorr(:,12),'--xk'  ) %plots for corrected
legend('Uncorrected','Corrected')
xlim([1500,5000])
title('T-motor Coplanar Torque Profile')
xlabel('RPM')
ylabel('Torque (Nm)')
hold off;

end