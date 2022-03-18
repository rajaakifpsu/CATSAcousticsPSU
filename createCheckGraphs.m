%Purpose: create graphs Greenwood asked for to check if the vold kelman 
%filter is working properly, ammong other things. Descriptions for each
%graph are commented in the code. Also calculates the PSD at each of the
%BPF harmonics for each mic.

%inputs:fpwelch is the vector of frequencies corresponding to other
%       variables.
%       originalSignal is the unfiltered pressure data. (nondeterministic
%       and deterministic.)
%       VKfiltered is the pressure data output the the vold kalman filter.
%       (deterministic)
%       residualBB is the pressure of broadband. (nondeterministic)
%       pxx is the pxx for unfiltered data. each row is a frequency, each
%       column is a mic.
%       PSD is 2d array representing power spectral density for unfiltered
%       pressure data for all the microphones. each row is a frequency,
%       eahc column a mic.
%       PSDBB is that but for broadband.
%       PSDFiltered is that but for the signal the vold kalman outputs
%       P_ref is reference pressure
%       PSDBBunsmoothed is like PSDBB but its calculated from the pxx for
%       broadband before it was run through the moving median filter.
%       PSDHPO is the PSD calculated from highpass only on the data
%       PSDBGnoPSU is psd of background noise without powere supply unit
%       PSDBGyesPSU is psd of background noise with power supply unit
%       PSDMN is psd of motor noise
%       fs is the sampling rate.
%       binwidth is width of window used in stuff. I dont think its
%       actually used in this code, though...
%       numblades is the number of blades
%       numofmics is the number of microphones
%       RPM is 2-element vector for RPM of rotor 1 and rotor 2.
%       savePath is the string or character array representing the folder
%       to save the plots to.
%       PSDBGnoPSU is the powe spectral density of background noise with no
%       power supply unit
%       PSDBGyesPSU is the same but with power supply unit
%       PSDMN is power spectral density of motor noise

%Outputs: PSDatBPFHarmonics is 2d array representing power spectral density 
%       at each of the BPF harmonics for each mic. each row is a mic, each
%       column is a harmonic of the BPF.
%       Also outputs a bunch of graphs. See the code comments for details.


%for testing: PSDatBPFHarmonics=createCheckGraphs(filteredTot, residualBBTot, pxxTot, PSDTot, PSDBBTot, P_ref, fs, binwidth, numblades, numofmics, RPM)

function PSDatBPFHarmonics=createCheckGraphs(fpwelch, originalSignal, VKfiltered, residualBB, pxx, PSD, PSDBB, PSDBBunsmoothed, PSDFiltered, PSDHPO, PSDBGnoPSU, PSDBGyesPSU, PSDMotorNoiseRed, P_ref, fs, binwidth, numblades, numofmics, RPM, savePath)

PSDMN=PSDMotorNoiseRed;

xlimits=[0, 20000]; %in hz

dt=1/fs; %time between samples
dfpwelch=fpwelch(2)-fpwelch(1);
times= dt* ( 0:(length(VKfiltered(:,1))-1) ); %times throughout whole test
BPFs=getBPFHarmonics(numblades,RPM(2));
freqIndexLims=getBPFfreqIndexLims(BPFs, fpwelch(2)-fpwelch(1), fpwelch, pxx(:,1)); %fpwelch(2)-fpwelch(1) is the difference in frequencies
angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert into radians (only for plotting)

x0=75; %distance from screen when graphing
y0=75; %same
width=1100; %graph size
height=400; %graph size
graphPos=[x0,y0,width,height];


%0) plots original signal v time, only for mic 1 (for simplicity)
figure()
plot(times,originalSignal)
title('Unfiltered signal v time (Mic 1)')
xlabel('time (s)')
ylabel('Pressure (Pa)')
saveDeletePlots(savePath);

for i=1:numofmics
    %1) plots Vkfiltered signal v time
    figure()
    plot(times,VKfiltered(:,i));
    title(sprintf('VK filtered signal v time (Mic %i)',i))
    xlabel('time (s)')
    ylabel('Pressure (Pa)')
    %1) plots PSDfiltered v fpwelch
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, sprintf('Vold-Kalman-Filtered PSD v frequency (Mic %i)',i), 1, xlimits, [0 100], graphPos, true)
end
saveDeletePlots(savePath);


%2) plots broadband v time, only for mic 1 again
figure()
plot(times,residualBB(:,1));
title('Broadband signal v time (Mic 1)')
xlabel('time (s)')
ylabel('Pressure (Pa)')
%2) plots all the signals v time
figure()
hold on;
plot(times, originalSignal(:,1))
plot(times, residualBB(:,1))
plot(times, VKfiltered(:,1))
title('Signals v time (Mic 1)')
xlabel('time (s)')
ylabel('Pressure (Pa)')
legend('Unfiltered Signal','Broadband Signal','VK filtered Signal')
hold off;
saveDeletePlots(savePath);
%Bellow graph is not useful. Done elsewhere.
% %2) plots broadband PSD v freq. note even though its using fpwelch,
% %thats fine. All the fpwelch vectors will be the same (assuming same
% %binwidth used, which it is... idk why i even made so many of them in main)
% createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'Broadband PSD v frequency', 1, xlimits, [0 100], graphPos, true)



%3) comparison of total, filtered, residual PSD (all mics)
for i=1:numofmics
    %PSD comparison H
    %first half of frequency range
    thistitle=['PSD Comparison H (Mic ',int2str(i),')'];
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, xlimits/2, [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, thistitle, i, xlimits/2, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    hold off;
    %now doing it for entire freq range
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    hold off;
    %now doing it for freq range of first 3 BPF harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, [0 BPFs(3)], [-100 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, [0 BPFs(3)], [-100 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, thistitle, i, [0 BPFs(3)], [-100 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    hold off;
    %now doing it for log scale on x axis.
    createPSDPlotLogScaleSingular(fpwelch, PSDBB, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSD, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    hold off;
    %Now doing it for x-axis scaled by the BPF harmonics, in first half of freq range
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, xlimits/2, [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, thistitle, i, xlimits/2, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    hold off;
    %Now doing the same but labelling them for harmonic number
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, xlimits/2, [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, thistitle, i, xlimits/2, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string((0:BPFs(1):(xlimits(2)/2)) / BPFs(1)))
    xlabel('BPF Harmonic')
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Broadband (without Peaks)','VK Filtered','Unfiltered','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD comparison A (BB with and without peaks. Motor noise too)
    %log scale
    thistitle=['PSD Comparison A (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDBB, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','Broadband (with Peaks)','Motor Noise','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDMN, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Broadband (without Peaks)','Broadband (with Peaks)','Motor Noise','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison B (BB w peaks and VK filter and motor noise)
    thistitle=['PSD Comparison B (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Broadband (with Peaks)','Vold-Kalman','Motor Noise','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDMN, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Broadband (with Peaks)','Vold-Kalman','Motor Noise','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison C (highpassonly and motor noise)
    thistitle=['PSD Comparison C (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDHPO, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Highpass Only','Motor Noise','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDHPO, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDMN, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Highpass Only','Motor Noise','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison D (background with psu, and motor noise)
    thistitle=['PSD Comparison D (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Background Noise with PSU','Motor Noise','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDMN, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Background Noise with PSU','Motor Noise','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison E (background with psu, background withou PSU)
    thistitle=['PSD Comparison E (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDBGnoPSU, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Background Noise with PSU','Background Noise without PSU','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBGnoPSU, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Background Noise with PSU','Background Noise without PSU','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison F (BB with peaks, Vk filter, highpass only, motor noise)
    thistitle=['PSD Comparison F (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDHPO, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Broadband (with Peaks)','Vold-Kalman','Highpass Only','Motor Noise','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'l', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDHPO, 'l', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDMN, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Broadband (with Peaks)','Vold-Kalman','Highpass Only','Motor Noise','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison G (altogether now)
    thistitle=['PSD Comparison G (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSD, 'this is filler', i, xlimits, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDHPO, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBB, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBGnoPSU, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, xlimits, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [0 60], graphPos, false)
    legend('Unfiltered','Highpass Only','Vold-Kalman','Broadband (with Peaks)','Broadband (without Peaks)','Background Noise without PSU','Background Noise with PSU','Motor Noise','Location','northeastoutside');
    hold off;
    %Now doing the same but zoomed in to first 10 harmonics
    createPowerSpectralDensityPlotSingular(fpwelch, PSD, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, true)
    hold on;
    createPowerSpectralDensityPlotSingular(fpwelch, PSDHPO, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDFiltered, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBGnoPSU, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, [0, BPFs(1)*10], [0 60], graphPos, false)
    createPowerSpectralDensityPlotSingular(fpwelch, PSDMN, thistitle, i, [0, BPFs(1)*10], [0 60], graphPos, false)
    legend('Unfiltered','Highpass Only','Vold-Kalman','Broadband (with Peaks)','Broadband (without Peaks)','Background Noise without PSU','Background Noise with PSU','Motor Noise','Location','northeastoutside');
    xticks(0:BPFs(1):(xlimits(2)/2))
    xticklabels(string(0:1:10))
    xlabel('BPF Harmonic')
    hold off;
    
    %PSD Comparison I (same as G, zoomed into a different y range)
    thistitle=['PSD Comparison I (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSD, 'this is filler', i, xlimits, [15 40], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDHPO, 'this is filler', i, xlimits, [15 40], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits, [15 40], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, xlimits, [15 40], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBB, 'this is filler', i, xlimits, [15 40], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBGnoPSU, 'this is filler', i, xlimits, [15 40], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, xlimits, [15 40], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits, [15 40], graphPos, false)
    legend('Unfiltered','Highpass Only','Vold-Kalman','Broadband (with Peaks)','Broadband (without Peaks)','Background Noise without PSU','Background Noise with PSU','Motor Noise','Location','northeastoutside');
    hold off;
    
    %PSD Comparison J (same as G, but zoomed into different x range)
    thistitle=['PSD Comparison J (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSD, 'this is filler', i, xlimits/2, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDHPO, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBBunsmoothed, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBB, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBGnoPSU, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDBGyesPSU, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits/2, [0 60], graphPos, false)
    legend('Unfiltered','Highpass Only','Vold-Kalman','Broadband (with Peaks)','Broadband (without Peaks)','Background Noise without PSU','Background Noise with PSU','Motor Noise','Location','northeastoutside');
    hold off;
    
    %PSD Comparison K. BB without peaks, VK filter, highpassonly, motor
    %noise
    thistitle=['PSD Comparison K (Mic ',int2str(i),')'];
    createPSDPlotLogScaleSingular(fpwelch, PSDBB, 'this is filler', i, xlimits/2, [0 60], graphPos, true)
    hold on;
    createPSDPlotLogScaleSingular(fpwelch, PSDFiltered, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDHPO, 'this is filler', i, xlimits/2, [0 60], graphPos, false)
    createPSDPlotLogScaleSingular(fpwelch, PSDMN, thistitle, i, xlimits/2, [0 60], graphPos, false)
    legend('Broadband (without Peaks)','Vold-Kalman','Highpass Only','Motor Noise','Location','northeastoutside');
    hold off;
    
    
    saveDeletePlots(savePath);
end



%4) looking at 1st, 11th and 12th mic to see if their peaks at BPF is the same.
createPowerSpectralDensityPlotSingular(fpwelch, PSD, 'Unfiltered PSD Near BPF', [1,11,12], [BPFs(1)-70 , BPFs(1)+70], [0 100], graphPos, true)
legend('mic 1','mic 11','mic 12');
% createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'Broadband PSD', [11,12], xlimits, [0 100], true)
% xlim([BPFs(1)-40 , BPFs(1)+40]);



%5) plot PSD zoomed in at lower range
createPowerSpectralDensityPlotSingular(fpwelch, PSD, 'Unfiltered PSD Before BPF', 1:(numofmics-1), [0 , BPFs(1)+50], [0 60], graphPos, true);
legend('Mic 1','Mic 2','Mic 3','Mic 4','Mic 5','Mic 6','Mic 7','Mic 8','Mic 9','Mic 10','Mic 11','Location','northeastoutside')
% createPowerSpectralDensityPlotSingular(fpwelch, PSDBB, 'Broadband PSD', 1:(numofmics-1), [0 , BPFs(1)+50], [0 60], graphPos, true);
% legend('Mic 1','Mic 2','Mic 3','Mic 4','Mic 5','Mic 6','Mic 7','Mic 8','Mic 9','Mic 10','Mic 11','Location','northeastoutside')



%6) Possible issue is that OASPL only slightly larger than SPL for BPF.
%Graphing SPL in range [0, BPF] should help identify the issue?
%Note thaat this will include the right side of the BPF tone. For none of
%that, see extra graph.
SPL=zeros(numofmics-1);
endindex=floor(BPFs(1)/dfpwelch)+1; %the index of frequency right before the BPF ( fpwelch(2) is the dfpwelch )
fprintf('BPF=%f     df=%f     endindex=%f',BPFs(1),fpwelch(2),endindex)
for i=1:(numofmics-1)
    SPL(i)=calculateSPL(fpwelch,pxx(:,i),P_ref,[1, endindex]);
end
figure()
polarplot(angle(1:(numofmics-1)),SPL(1:(numofmics-1)),'k')
thetalim([-40,30]);
rlim([40,80]);
title('Unfiltered SPL [0-BPF]');
saveDeletePlots(savePath);



%extra) plot for SPL of [0, left limit of BPF tone], and [Right limit of
%BPF tone, 20000], and [0, 20000 (aka OASPL)]
SPLbeforeBPFTone=zeros(1,(numofmics-1)); %initializing
SPLafterBPFTone=SPLbeforeBPFTone;
OASPL=SPLbeforeBPFTone;
for i=1:(numofmics-1)
    SPLbeforeBPFTone(i)=calculateSPL(fpwelch,pxx(:,i),P_ref,[1,freqIndexLims(1,1)]);
    SPLafterBPFTone(i)=calculateSPL(fpwelch,pxx(:,i),P_ref,[freqIndexLims(2,1),length(pxx(:,1))]);
    OASPL(i)=calculateSPL(fpwelch,pxx(:,i),P_ref,[1,length(pxx(:,1))]); %could probably just pass this into this function from elsewhere, but this may be good check if code working or not
end
figure()
polarplot(angle(1:(numofmics-1)),SPLbeforeBPFTone,'g',angle(1:(numofmics-1)),SPLafterBPFTone,'r',angle(1:(numofmics-1)),OASPL,'k')
thetalim([-40,30]);
rlim([30,80]);
title('Unfiltered SPL of Various Frequency Ranges');
legend('Before BPF Tone','After BPF Tone','OASPL')



%table) a table for the values of PSD and the BPF and its harmonics. 
PSDatBPFHarmonics=zeros(numofmics,length(BPFs)); %each column a harmonic, each row a mic
for i=1:numofmics %loops through all mics
    for j=1:length(BPFs) %loops through all harmonics of bpf
        %accounts for difference between index of frequency vector before harmonic and the frequency of the actual harmonic
        BPFindex=floor(BPFs(j)/dfpwelch)+1; %the index of frequency right before the harmonic
        slope=(PSD(BPFindex+1,i)-PSD(BPFindex,i))/dfpwelch; %slope in the PSD-frequency graph
        PSDatBPFHarmonics(i,j)=PSD(BPFindex,i)+((BPFs(j)-(fpwelch(BPFindex)))*slope);
    end
end


end