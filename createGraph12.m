%fig is a 2x3 array for which figures to use. 6 total.
%groupLength is the number of tests in a group. (ie in the vdhcsab group as opposed to the vdhcsbb)

function createGraph12(avFMRPMCorr,SEMCorr, avFMRPMCorrTOT, SEMCorrTOT, groupLength,fig)
loop=1;
for i=1:groupLength:length(avFMRPMCorr(:,1)) %loops through each group of 8 tests. Ie the VDHCSAB vs VDHCSBB vs VDHCSCB groups.
    figure(fig(loop,1))
    grid on;
    hold on;
    errorbar(avFMRPMCorr(i:i+groupLength-1,13),avFMRPMCorr(i:i+groupLength-1,3),SEMCorr(i:i+groupLength-1,3),'--xb'  ) %plots Fz v RPM for corrected motor1
    errorbar(avFMRPMCorr(i:i+groupLength-1,14),avFMRPMCorr(i:i+groupLength-1,9),SEMCorr(i:i+groupLength-1,9),'--xr'  ) %same for motor2
    errorbar(avFMRPMCorrTOT(i:i+groupLength-1,7),avFMRPMCorrTOT(i:i+groupLength-1,3),SEMCorrTOT(i:i+groupLength-1,3),'--xk') %same for total
    legend('Motor 1','Motor 2','Total','Location','northwest')
    xlim([1500,5100])
    ylim([0,60])
    title('Seperated Harmony Straight Thrust Profile')
    xlabel('RPM')
    ylabel('Thrust (N)')
    hold off;

    figure (fig(loop,2))
    grid on;
    hold on;
    errorbar(avFMRPMCorr(i:i+groupLength-1,14),avFMRPMCorr(i:i+groupLength-1,6),SEMCorr(i:i+groupLength-1,6),'--xb'  ) %plots Mz v RPM for corrected motor1
    errorbar(avFMRPMCorr(i:i+groupLength-1,14),avFMRPMCorr(i:i+groupLength-1,12),SEMCorr(i:i+groupLength-1,12),'--xr'  ) %same for motor2
    errorbar(avFMRPMCorrTOT(i:i+groupLength-1,7),avFMRPMCorrTOT(i:i+groupLength-1,6),SEMCorrTOT(i:i+groupLength-1,6),'--xk') %same for total
    legend('Motor 1','Motor 2','Total','Location','southwest')
    xlim([1500,5100])
    ylim([-2.25,0])
    title('Seperated Harmony Straight Torque Profile')
    xlabel('RPM')
    ylabel('Torque (Nm)')
    hold off;
    
    loop=loop+1;
end
end