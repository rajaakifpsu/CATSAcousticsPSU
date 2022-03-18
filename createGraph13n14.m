%figs is 2x3 matrix for which figures to use
%ylimit is a 3x4 matrix for what ylimits to use. first 2 columns are for force, 2nd two are for moment.
%groupLen is number of tests in a group. ie the vdhcsab vs vdhcsbb groups

function createGraph13n14(avFMRPMCorr,SEMCorr, avFMRPMCorrTOT, SEMCorrTOT, groupLen, fig, ylimit)
%title('T-motor Coplanar Force Profiles') <------Add title in powerpoint
outerloop=1;
for j=1:groupLen:24 %loops through each group of 8 tests. Ie the VDHCSAB vs VDHCSBB vs VDHCSCB groups.
    axistitles=["Fx (N)", "Fy (N)","Fz (N)"];
    figure(fig(outerloop,1))
    for i=1:3 %loop through Fx, Fy, Fz
        subplot(3,1,i)
        hold on;
        grid on;
        errorbar(avFMRPMCorr(j:j+groupLen-1,13),avFMRPMCorr(j:j+groupLen-1,i),SEMCorr(j:j+groupLen-1,i),'--xb'  ) %plots for corrected in motor1
        errorbar(avFMRPMCorr(j:j+groupLen-1,14),avFMRPMCorr(j:j+groupLen-1,i+6),SEMCorr(j:j+groupLen-1,i+6),'--xr'  ) %same for motor2
        errorbar(avFMRPMCorrTOT(j:j+groupLen-1,7) , avFMRPMCorrTOT(j:j+groupLen-1,i),SEMCorrTOT(j:j+groupLen-1,i),'--xk' ) %same for total
        legend('Motor 1','Motor 2', 'Total','Location','northeastoutside')
        xlim([2500,5100])
        ylim(ylimit(outerloop,1:2))
        xlabel('RPM')
        ylabel(axistitles(i))
        hold off;
    end

    figure(fig(outerloop,2))
    %title('T-motor Coplanar Moment Profiles') <-----Add title in powerpoint
    axistitles=["Mx (Nm)","My (Nm)","Mz (Nm)"];
    for i=1:3 %loop through Mx, My, Mz
        subplot(3,1,i)
        hold on;
        grid on;
        errorbar(avFMRPMCorr(j:j+groupLen-1,13),avFMRPMCorr(j:j+groupLen-1,i+3),SEMCorr(j:j+groupLen-1,i+3),'--xb'  ) %plots corracted in motor1
        errorbar(avFMRPMCorr(j:j+groupLen-1,14),avFMRPMCorr(j:j+groupLen-1,i+9),SEMCorr(j:j+groupLen-1,i+9),'--xr'  ) %same for motor2
        errorbar(avFMRPMCorrTOT(j:j+groupLen-1,7) , avFMRPMCorrTOT(j:j+groupLen-1,i+3),SEMCorrTOT(j:j+groupLen-1,i+3),'--xk' ) %same for total
        legend('Motor 1','Motor 2','Total','Location','northeastoutside')
        xlim([2500,5100])
        ylim(ylimit(outerloop,3:4))
        xlabel('RPM')
        ylabel(axistitles(i))
        hold off;
    end
    outerloop=outerloop+1;
end

end