function createGraph6(avFMRPM,targettedRPMs, fig)

figure(fig)
hold on;
grid on;
scatter(targettedRPMs,(avFMRPM(:,14)'-targettedRPMs),40,'k','x')
xlabel('Target RPM')
ylabel('Measured RPM - Target RPM')
xlim([1400,5000])
ylim([-25,25])
title('RPM Error')
hold off;

end