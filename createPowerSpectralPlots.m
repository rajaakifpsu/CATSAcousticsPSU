function createPowerSpectralPlots(fpwelch, PSD, graphtitles)

figure()
plot(fpwelch,PSD)
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density -Full Data (dB)')
xlim([0 20000])
ylim([0 100])
title(graphtitles(1))
legend('mic 1', 'mic 2','mic 3','mic 4','mic 5','mic 6','mic 7','mic 8','mic 9','mic 10','mic 11','mic 12','Location','northeastoutside')

%Does the same for first and last mic
figure()
hold on;
plot(fpwelch,PSD(:,1));
plot(fpwelch,PSD(:,length(PSD(1,:))));
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB)')
xlim([0 20000])
ylim([0 100])
title(graphtitles(2))
legend('mic 1', 'mic 12','Location','northeastoutside')
hold off;


for i=1:numofmics %does the same but for each mic individually
    figure()
    plot(fpwelch,PSD(:,i))
    xlabel('Frequency (Hz)')
    ylabel('Power Spectrum (dB)')
    xlim([0 20000])
    ylim([0 100])
    title([graphtitles(3),int2str(i)])
end

end