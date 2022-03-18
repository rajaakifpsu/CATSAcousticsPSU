%Purpose: to create power spectral density plots in multiple configurations

%inputs: fpwelch is vector of frequencies corresponding to PSD
%       PSD is a matrix representing the power spectral density of the
%       test. each row corresponds to a frequency, each column to a mic.
%       can be the broadband, the vk, whatever the function is passed.
%       numofmics is the numbeer of microphones
%       numblades is the number of blades
%       RPMmotor2 is the average RPM of the second motor (the one rotating)
%       graphtitles is a string array for the titles of the graphs to be
%       made. 
%       savePath is the string or character array for the folder to save
%       plots into.

%outputs: a bunch of plots, but no new data.

function createPowerSpectralDensityPlots(fpwelch, PSD, numofmics, numblades, RPMmotor2, graphtitles, savePath)

xlimits=[0 20000]; % in Hz
ylimits=[0 100];

BPFs=getBPFHarmonics(numblades,RPMmotor2); %the blade pass frequency

x0=75; %distance from screen when graphing?
y0=75; %same
width=1100; %graph size
height=400; %graph size
graphPos=[x0,y0,width,height];

%plots alll mics together
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(1), 1:numofmics, xlimits, ylimits, graphPos, true)
legend('mic 1', 'mic 2','mic 3','mic 4','mic 5','mic 6','mic 7','mic 8','mic 9','mic 10','mic 11','mic 12','Location','northeastoutside')
%does all mics together but with log scaling
createPSDPlotLogScaleSingular(fpwelch, PSD, graphtitles(1), 1:numofmics, xlimits, ylimits, graphPos, true)
legend('mic 1', 'mic 2','mic 3','mic 4','mic 5','mic 6','mic 7','mic 8','mic 9','mic 10','mic 11','mic 12','Location','northeastoutside')
%does all mics together but with different scaling (scaled by the BPF), and only for half the span of previous
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(1), 1:numofmics, xlimits/2, ylimits, graphPos, true)
legend('mic 1', 'mic 2','mic 3','mic 4','mic 5','mic 6','mic 7','mic 8','mic 9','mic 10','mic 11','mic 12','Location','northeastoutside')
xticks(0:BPFs(1):xlimits(2)/2)
%does the same but labels the ticks by the number of the harmonic
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(1), 1:numofmics, xlimits/2, ylimits, graphPos, true)
legend('mic 1', 'mic 2','mic 3','mic 4','mic 5','mic 6','mic 7','mic 8','mic 9','mic 10','mic 11','mic 12','Location','northeastoutside')
xticks(0:BPFs(1):xlimits(2)/2)
xticklabels(string( (0:BPFs(1):xlimits(2)/2)/BPFs(1) ))
xlabel('BPF Harmonic')
%does the same but zoomed in on first 10 harmonics
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(1), 1:numofmics, [0, BPFs(1)*10], ylimits, graphPos, true)
legend('mic 1', 'mic 2','mic 3','mic 4','mic 5','mic 6','mic 7','mic 8','mic 9','mic 10','mic 11','mic 12','Location','northeastoutside')
xticks(0:BPFs(1):(BPFs(1)*10))
xticklabels(string( (0:1:10 )))
xlabel('BPF Harmonic')

%Does the first and last mic only
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(2), [1,numofmics], xlimits, ylimits, graphPos, true)
legend('mic 1', 'mic 12','Location','northeastoutside')
%Does the first and last mic but with log scaling
createPSDPlotLogScaleSingular(fpwelch, PSD, graphtitles(2), [1,numofmics], xlimits, ylimits, graphPos, true)
legend('mic 1', 'mic 12','Location','northeastoutside')
%does the same but with different scaling (scaled by the BPF), and only for half the span of previous
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(2), [1,numofmics], xlimits/2, ylimits, graphPos, true)
legend('mic 1', 'mic 12','Location','northeastoutside')
xticks(0:BPFs(1):xlimits(2)/2)
%does the same but labels the ticks by number of harmonic
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(2), [1,numofmics], xlimits/2, ylimits, graphPos, true)
legend('mic 1', 'mic 12','Location','northeastoutside')
xticks(0:BPFs(1):xlimits(2)/2)
xticklabels(string( (0:BPFs(1):xlimits(2)/2)/BPFs(1) ))
xlabel('BPF Harmonic')
%does the same but zoomed in on first 10 harmonics
createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitles(2), [1,numofmics], [0, BPFs(1)*10], ylimits, graphPos, true)
legend('mic 1', 'mic 12','Location','northeastoutside')
xticks(0:BPFs(1):(BPFs(1)*10))
xticklabels(string( 0:1:10 ))
xlabel('BPF Harmonic')

saveDeletePlots(savePath);

for i=1:numofmics %does each mic individually
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, [graphtitles(3),int2str(i)], i, xlimits, ylimits, graphPos, true) %normal scaling
    
    createPSDPlotLogScaleSingular(fpwelch, PSD, [graphtitles(3),int2str(i)], i, xlimits, ylimits, graphPos, true) %log scaling
    
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, [graphtitles(3),int2str(i)], i, xlimits/2, ylimits, graphPos, true)
    xticks(0:BPFs(1):xlimits(2)/2) %different scaling, by BPF
    
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, [graphtitles(3),int2str(i)], i, xlimits/2, ylimits, graphPos, true)
    xticks(0:BPFs(1):xlimits(2)/2)
    xticklabels(string((0:BPFs(1):xlimits(2)/2) / BPFs(1))) %different scaling, by BPF Harmonic
    xlabel('BPF Harmonic')
    
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, [graphtitles(3),int2str(i)], i, [0, BPFs(1)*10], ylimits, graphPos, true)
    xticks(0:BPFs(1):(BPFs(1)*10))
    xticklabels(string(0:1:10)) %different scaling, by BPF Harmonic, zoomed in on first 10 harmonics
    xlabel('BPF Harmonic')
end


end