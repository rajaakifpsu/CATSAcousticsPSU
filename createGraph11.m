%fig is a 2x3 array for which figures to use. 6 total.


function createGraph11(avFMRPMCorr, avFMRPMCorrTOT, targettedRPMs, groupLen, fig)
loop=1;
for i=1:groupLen:24 %loops through each group of 8 tests. Ie the VDHCSAB vs VDHCSBB vs VDHCSCB groups.
    figure(fig(loop))
    hold on;
    grid on;
%     disp(avFMRPMCorr(i:i+7,13))
%     disp(targettedRPMs(i:i+7))
    scatter(targettedRPMs(i:i+groupLen-1),(avFMRPMCorr(i:i+groupLen-1,13)'-targettedRPMs(i:i+groupLen-1)),40,'b','x')
    scatter(targettedRPMs(i:i+groupLen-1),(avFMRPMCorr(i:i+groupLen-1,14)'-targettedRPMs(i:i+groupLen-1)),40,'r','x')
    scatter(targettedRPMs(i:i+groupLen-1),(avFMRPMCorrTOT(i:i+groupLen-1,7)'-targettedRPMs(i:i+groupLen-1)),40,'k','x')
    xlabel('Target RPM')
    ylabel('Measured RPM - Target RPM')
    xlim([1400,5100])
    ylim([-35,35])
    title('RPM Error')
    legend('Motor 1','Motor 2','Total','Location','northeastoutside')
    hold off;
    
    loop=loop+1;
end
end