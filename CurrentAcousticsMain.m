%TO DO: still have roll off in the residual. Need to fix this somehow.

clear all;
close all ;
clc;

tic

%
numofmics = 12; 

% % % % %USE IF YOU HAVE THE .MAT FILE
% % % % 
% % % % % load('tdhcs1.mat');
% % % % % numofmicss = 1;
% % % % % SensitivityData = SensitivityData./1000; %convert mV to V
% % % % % for j = 1:12
% % % % %     Pressure(:,j) = AcousticData(:,j)/SensitivityData(j);
% % % % % end
% % % % % N = length(Pressure); %length of signal
% % % % % dt = 1/fs;    %time between each input
% % % % % time = dt*[0:N-1];  %get array of times at each point
% % % % % T = max(time); %max time or period
% % % % % df = 1/T; 

% % % % %use if you have BK file
% % % % %filename = 'AnechoicTest_kdhcs-9.h5';
% % % % %data = readBKHDF5(filename); %Read BK file
% % % % %fs_total = data.sampling_rate; %sampling rate of all mics
% % % % %max_time = max(data.Signals{1,1}(:,1)); %length of acoustic signal (s)
% % % % %time_array = (data.Signals{1,1}(:,1)); %time array of the run
% % % % %T = max_time;

%Use if you dont have the .mat file & only have the Labview System H5

%%%sepdata = readh5('acs_data_tdhcs7.h5'); %reading the HDF5 file (acoustics)
%%%sepdata1 = readmatrix('tdhcs7.txt'); %performance data

filename = "acs_data_tdhcs7.h5";
file = hdf5info(filename);
perfFileName="tdhcs7.txt"; %performance
performanceMatrix=readmatrix(perfFileName);

%Reading all the Signals from the h5 file and calibrating the signals to Pa
%to work with the processing techniques
fs = hdf5read(file.GroupHierarchy.Datasets(4)); %sampling rate in Hz. Which is the 4th field
acousticsdata_uncorr = hdf5read(file.GroupHierarchy.Datasets(1)); % Acoustics value in voltage.
micsensitivities = hdf5read(file.GroupHierarchy.Datasets(5)); % Microphone sensitivities mv/Pa
acousticsdata = (1000./(micsensitivities')).*acousticsdata_uncorr; % Calibrated Acoustic values%initialize calculations
noofsamples = length(acousticsdata); %length of signal

N = noofsamples;
dt = 1/fs;    %time increment between each data point
time_array = [0:N-1]*dt;  %get array of times at each point
T = max(time_array); %max time or period
df = 1/T; 

%% take all the pressure time data here
for i = 1:numofmics
    acoustics_total(:,i) = acousticsdata(:,i); %raw sound pressures
    pressure_total(:,i) = acoustics_total(:,i); 
end 


%%  taking a section of time in the .h5 file

starttime = 90; %start time in seconds
endtime = 120; %end time in seconds
T_reduced = endtime - starttime;
N_reduced = T_reduced*fs;

for i = 1:numofmics %for every microphone
    
        acoustics_reduced(:,i) = acousticsdata((starttime*fs:endtime*fs-1),i); %sound pressures
        time_array_reduced(:,i) = time_array(starttime*fs:endtime*fs-1); %getting the reduced time for each mic
        pressure_reduced(:,i) = acoustics_reduced(:,i);
end

%% Correct the DC Offset

for i = 1:numofmics
   
   meanDCoffset(:,i) = mean(pressure_total(:,i),1);
   pressure_total(:,i) = pressure_total(:,i) - meanDCoffset(1,i);
   pressure_reduced(:,i) = pressure_reduced(:,i) - meanDCoffset(1,i);
    
end

%% calculate the OASPL RMS

 dbref = 2e-5; % reference pressure in Pa or 2e-5;
 P_ref = dbref;
 G_ref = dbref^2; 

for k = 1:numofmics
    rmspress_total(k,:) = rms(pressure_total(1:N,k));
    OASPL_total_RMS(k,:) = 20*log10(rmspress_total(k,:)/P_ref); %OASPL at each Mic
    OASPL_reduced_RMS(k,:) = 20*log10(rms(pressure_reduced(1:N_reduced,k))/P_ref); %OASPL reduced at each Mic
    deltatotalreduced = OASPL_total_RMS - OASPL_reduced_RMS;
end 


%% Angles of the Microphone Array

angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert angles to radians

horzdistance_centerMic = 60*0.0254; %inches to meters
vertdistance_centerMic = 78*0.0254; %inches to meters

%% define RPM
performanceSampRate=20000;
RPM = getRPMSteady(performanceMatrix, performanceSampRate);
numblades = 2;
%% P welch method
 bandwidth = 20; % width of the record in frequency (Hz)
 binwidth = 2^nextpow2(fs/bandwidth); %not necesary for math, but reduces time needed to run code.
 
 %doing pwelch on the full set of data
% % % % % [pxxTot, fpwelchTot, PSDTot,PowerSpectrumTot, energyTot, energytotTot, oasplTot, pxxBBTot, fpwelchBBTot, PSDBBTot, PowerSpectrumBBTot, energyBBTot, enrgytotBBTot, oasplBBTot, filteredTot, residualBBTot] = pwelchFullAndBB(pressure_total, binwidth,fs,P_ref,numofmics, RPM, numblades); %OLD SHIT: pwelchFullAndBB(pressure_total, binwidth,fs,P_ref,numofmics);   
%now on the reduced set of data
[pxxRed, fpwelchRed, PSDRed,PowerSpectrumRed, energyRed, energytotRed, oasplRed, pxxBBRed, fpwelchBBRed, PSDBBRed, PowerSpectrumBBRed, energyBBRed, enrgytotBBRed, oasplBBRed, filteredRed, residualBBRed] = pwelchFullAndBB(pressure_reduced, binwidth,fs,P_ref,numofmics, RPM, numblades);

%% Plots
% % % % % %Graphs power spectrum aka narrowband spectrum vs frequency (all 12 mics)
% % % % % graphTitlesTot=["Total PSD for the whole range of points","Total PSD for the whole range of points","Total PSD for microphone "];
% % % % % createPowerSpectralDensityPlots(fpwelchTot, PSDTot, numofmics, graphTitlesTot);
% % % % % 
% % % % % graphTitlesBB=["Broadband PSD for the whole range of points","Broadband PSD for the whole range of points","Broadband PSD for microphone "];
% % % % % createPowerSpectralDensityPlots(fpwelchBBTot, PSDBBTot, numofmics, graphTitlesBB);
% % % % % 
% % % % % %Directivity
% % % % % 
% % % % % angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
% % % % % angle = angle.*(pi/180); %convert into radians (only for plotting)
% % % % % createDirectivityPlot(oasplTot,oasplBBTot,numofmics,'Full-Data Directivity Plot')
% % % % % createDirectivityPlot(oasplRed,oasplBBRed,numofmics,'Reduced-Data Directivity Plot')
% % % % % figure()
% % % % % polarplot(angle(1:(numofmics-1)),oasplBBTot(1:(numofmics-1)),'r',angle(1:(numofmics-1)),oasplBBRed(1:(numofmics-1)),'k') %the broadband for total and reduced data
% % % % % thetalim([-40 , 30]);
% % % % % rlim([60 , 80]);
% % % % % legend('Total Data','Reduced Data')
% % % % % title('Broadband Directivity Plot')
% % % % % 
% % % % % %Directivity of BPF and Harmonics therof
% % % % % SPLForBPFHarmonicsTot=createDirectivityPlotForHarmonics(pxxTot,oasplTot,fpwelchTot,RPM(2),numblades,numofmics,P_ref,'BPF Harmonics SPL Directivity Plot (Total Data)');
% % % % % SPLForBPFHarmonicsRed=createDirectivityPlotForHarmonics(pxxRed,oasplRed,fpwelchRed,RPM(2),numblades,numofmics,P_ref,'BPF Harmonics SPL Directivity Plot (Reduced Data)');
% % % % % 
% % % % % 
% % % % % %Graphs Greenwood asked for to see if VK works correctly
% % % % % createCheckGraphs(filteredTot, residualBBTot, pxxTot, PSDTot, PSDBBTot, P_ref, fs, binwidth, numblades, numofmics, RPM)
% % % % % 
% % % % % %Graphs for NASA
% % % % % createNasaDeliverables(fpwelchTot, pxxTot,SPLForBPFHarmonics, P_ref, angle, numofmics);

%FOR TESTING purposes: comment out later...
nooverlap=binwidth*0.5; %50 percent overlap
nfft=binwidth; %number of points to use in the fft=binwidth.
dfForPwelch=fs/binwidth; 
[W, ~, ~] = Hann_window(binwidth); %creates hann W sized to binwidth
ef = sum(W.*W)/binwidth; % mean-square of W & 1/Aw = energy_factor 

[pxxFilteredRed ,fpwelchFiltered]= pwelch(filteredRed,W,nooverlap,[],fs, 'onesided');
PSDFilteredRed = 10*log10(pxxFilteredRed/(P_ref^2)); %psd (Pa^2/Hz) 
createPowerSpectralDensityPlotSingular(fpwelchRed, PSDRed, 'this is filler', 1, [0 20000], [0 100], true)
hold on;
createPowerSpectralDensityPlotSingular(fpwelchRed, PSDFilteredRed, 'this is filler', 1, [0 20000], [0 100], false)
createPowerSpectralDensityPlotSingular(fpwelchRed, PSDBBRed, 'Reduced Data PSD Comparison', 1, [0 20000], [0 100], false)
legend('Unfiltered','VK Filtered','Broadband','Location','northeastoutside');
hold off;


toc






% % % 
% % % for k = 1:numofmics
% % % 
% % %     fq = 0:bandwidth:fs/2;
% % %     %[pxx,f] = pwelch(pressure_total(:,k) ,W,nooverlap,[],fs);
% % %     %[pxx(:,k),f(:,k)] = pwelch(pressure_total(:,k),W,nooverlap,[],fs);
% % %     %W is a vector
% % %     [pxx(:,k),fq(:,k)] = pwelch(pressure_total(:,k),80000,512,[],fs); %80000 is interget (where issue occurs)
% % %     %[pxx(:,k),f(:,k)] = pwelch(pressure_reduced(:,k),W,nooverlap,[],fs);
% % %     
% % %     
% % %     %[pxx,f] = pwelch(x,W,nooverlap,f,fs) returns the two-sided Welch PSD estimates at the frequencies specified in the vector, f. The vector f must contain at least two elements, because otherwise the function interprets it as nfft. The frequencies in f are in cycles per unit time. The sample rate, fs, is the number of samples per unit time. If the unit of time is seconds, then f is in cycles/sec (Hz).
% % %     
% % %     %from miccal
% % % %     Spl_denorm = pxx / bandwidth; 
% % % %     energy_tot = sum(Spl_denorm);
% % % %     oaspl = 10*log10((energy_tot)/Pref.^2);
% % % 
% % %     
% % %     %SPL(:,i) = 10*log10(Gxx(:,i)./Gref);
% % %     %SPLp(:,k) = 10*log10(pxx(:,k)./G_ref); %SPL at each Mic
% % %     SPLp(:,k) = 10*log10(pxx(:,k)./G_ref/bandwidth); %SPL at each Mic
% % % 
% % %     %valp(:,k) = (1/T)*(sum((pxx(:,k)))); %energy in the system
% % %     %OASPLp(:,k) = 10*log10(valp(:,k)./G_ref); %OASPL at each Mic  
% % %    
% % %     
% % %     Spl_denorm(:,k) = pxx(:,k) / bandwidth; 
% % %     energy_tot(:,k) = sum(Spl_denorm(:,k));
% % %     oaspl(:,k) = 10*log10((energy_tot(:,k))/P_ref.^2);
% % %     %RMS_Gxx(:,k) = sqrt(sum((Gxx_use(:,k)))*df); %root mean squared energy. or root of energy.
% % %     %OASPL_exp_gxxmethod(:,k) = 20*log10(RMS_Gxx(:,k)./P_ref);
% % %     
% % %     RMS_pxx(k) = sqrt(sum(pxx(:,k))*df); %root mean squared energy. or root of energy.
% % %     OASPL_exp_pxxmethod(k) = 10*log10(RMS_pxx(k)./G_ref);
% % % 
% % % 
% % %     figure()
% % %     hh = plot(f(1:nooverlap-1,k), pxx(1:nooverlap-1,k), '-b'); %skip point 1
% % %     set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% % %     set(gca, 'fontsize', 13) % gca = "get current axis"
% % %     grid on
% % %     xlim([0 20000]) % x limit to half of samplerate
% % %     title({'Gxx for mic ',num2str(k)})
% % %     xlabel('Frequency [Hz]','fontsize', 13)
% % %     ylabel('Gxx [Pa^{2}/Hz]', 'fontsize', 13)
% % % 
% % %     figure()
% % %     hh = plot(f(61:nooverlap-1,k), SPLp(61:nooverlap-1,k), '-b'); %skip point 1
% % %     set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% % %     set(gca, 'fontsize', 13) % gca = "get current axis"
% % %     grid on
% % %     xlim([0 20000]) % x limit to half of samplerate
% % %     title({'SPL for mic ',num2str(k)})
% % %     xlabel('Frequency [Hz]','fontsize', 13)
% % %     ylabel('SPL [dB]', 'fontsize', 13)
% % % 
% % % end
% % %     
% % %     angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
% % %     angle = angle.*(pi/180);
% % %     figure()
% % %     rlim = [30 45]; %dBs
% % %     tlim = [-36.62 23.01];  %degrees
% % %     polarplot(angle(1:11),OASPLp(1:11),'LineWidth',1.3)
% % %     ax=gca;
% % %     ax.FontSize=14;
% % %     grid on
% % %     grid minor
% % %     ax.RLim=rlim;
% % %     ax.ThetaLim=tlim;
% % %     legend({'w/Filter','Original'})
% % %     title(['OASPL Directivity for ',num2str(RPM),...
% % %           'rpm [Hover]'])
% % %     % saveas(gca,['BLtd Plot for Sep ',num2str(RPM),'rpm[Edgewise].png'])
% % %     % saveas(gca,['Tonal Directivity for Sep ',num2str(RPM),'rpm[Edgewise].png'])
% % % 
% % %     
%% For each Mic to check the value and learn without using any filter or overlaps
% 
% freq_axis = ((1:N)*df)';
%  
% lsp = fft(acoustics_total)*dt; %Find the linear spectrum, lsp 
% Sxx = abs(lsp).^2./T; %get dual sided spectrum,[Pa^2/Hz]. Here Sxx size is numofsamples
% 
%set requirements
% %at N/2| Gxx(N/2) = Sxx(N/2)
% %at 0| Gxx(1) = Sxx(1)
% %others| Gxx = 2*Sxx from 1 to N/2 - 1 in freq domain
% 
% %for mic 1
% mic = 1;
% Sxx_mic = Sxx(:,mic);
% Gxx = 2*Sxx_mic;  %(2/T)*abs(lsp).^2
% Gxx(1) = Sxx_mic(1);
% Gxx(N/2+1) = Sxx_mic(N/2+1); 
% %Gxx has Noofsamples size
% %Sxx has Noofsamples size
% 
% % Compare Mean of Time Series with the convention
% % conventional must = to Ms time series
% 
% % time domain version. Mean-squared value of the original time series.
% MS_timeSeries = sum(acoustics_total(:,1).^2)./N;
% 
% % Sxx*df
% MS_Xm = sum(Sxx_mic*df); % integral of double sided. conventional
% 
% % Gxx*df but Gxx is to half of number of samples
% MS_Gxx = sum(Gxx(1:(N/2+1))*df); %integral of single sided. conventional
% 
% % if the MS time series =/= Sxx*df to all number of samples, Sxx is wrong
% % if the MS time series =/= Gxx*df to half number of samples, Gxx is wrong
% 
% % plot the Gxx to half of the number of samples
% figure()
% hh = plot(freq_axis(2:(N/2)+1), Gxx(2:(N/2)+1), '-b'); %skip point 1
% set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% set(gca, 'fontsize', 13) % gca = "get current axis"
% grid on
% xlim([0 40000]) % x limit to half of samplerate
% title({'Gxx for mic ',num2str(mic)})
% xlabel('Frequency [Hz]','fontsize', 13)
% ylabel('Gxx [Pa^{2}/Hz]', 'fontsize', 13)
% 
% %get log plot
% figure()
% hh = semilogy(freq_axis(2:(N/2)), Gxx(2:(N/2)+1), '-b'); %skip point 1
% set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% set(gca, 'fontsize', 13) % gca = "get current axis"
% xlim([0 40000]) % x limit to half of samplerate
% grid on
% title({'Log plot of Gxx for mic ',num2str(mic)})
% xlabel('Frequency [Hz]','fontsize', 13)
% ylabel('Gxx [Pa^{2}/Hz]', 'fontsize', 13)

%%
