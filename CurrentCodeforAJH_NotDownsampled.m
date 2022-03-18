% %TO DO MANUALLY: will have to input the temp and pressure vectors. can
% %probably make this all into a function that takes these as perameters.
% 
% %TO CHANGE WHEN RUNNING VS WHEN TESTING: numofmics from =12 to =1.
% 
clear all;
close all ;
clc;

tic

%flyover condition 
%assuming all the microphones have the same sampling rate

%% Loading in data for performance and acoustics

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

filename = "acs_data_tdhcs7.h5"; %acoustics
file = hdf5info(filename);
perfFileName="tdhcs7.txt"; %performance
performanceMatrix=readmatrix(perfFileName);
load('NoisePSDs.mat') %the background noise (with and without psu) and the
%motor noise. They are calculated elsewhere in another code,
%CurrentCodeforBackgroundNoise.m

%Reading all the Signals from the h5 file
%to work with the processing techniques
fs = hdf5read(file.GroupHierarchy.Datasets(4)); %sampling rate in Hz. Which is the 4th field
acousticsdata_uncorr = hdf5read(file.GroupHierarchy.Datasets(1)); % Acoustics value in voltage.

savePath = 'C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_17_2022_AJH';
% the location to save figures to.

%% NOT DOING THIS
% % % %% Downsampling the signals
% % % 
% % % acousticsdata_uncorr = acousticsdata_uncorr(1:2:length(acousticsdata_uncorr(:,1)),:); %downsampling the signals for computational efficiency
% % % fs=fs/2; %downsampling the signal for computational efficiency

%% Converting from voltage to Pa with microphone sensitivities

micsensitivities = hdf5read(file.GroupHierarchy.Datasets(5)); % Microphone sensitivities mv/Pa
acousticsdata = (1000./(micsensitivities')).*acousticsdata_uncorr; % Calibrated Acoustic values%initialize calculations
noofsamples = length(acousticsdata); %length of signal

%% Correct the distance from rotor to microphones to be the distance NASA wants.

correctDist=1.8956; %in m (the distance NASA wants)
micDist=[65.19, 62.97, 61.34, 60.34, 60.00, 60.34, 61.34, 62.97, 65.19, 67.93, 71.14, 74.75];%in inches. distance we have to each mic in test.
micDist=convlength(micDist, 'in','m'); %converts micDist to m

for i=1:numofmics
    acousticsdata(:,i) = acousticsdata(:,i) * (micDist(i)/correctDist); %its just signal multiplied by (dist we have / dist we want)
end

%% Defining other little variables

N = noofsamples;
dt = 1/fs;    %time increment between each data point
time_array = [0:N-1]*dt;  %get array of times at each point
T = max(time_array); %max time or period
df = 1/T; 

%% take all the pressure time data here

for i = 1:numofmics
    acoustics_total(:,i) = acousticsdata(:,i); 
    pressure_total(:,i) = acoustics_total(:,i); 
end 

%%  taking a section of time in the .h5 file to define pressure_reduced

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
%%%[pxxTot,fpwelchTot, PSDTot,PowerSpectrumTot, oasplTot, pxxBBTot, PSDBBTot, PowerSpectrumBBTot, oasplBBTot, pxxVKTot, PSDVKTot, PowerSpectrumVKTot, oasplVKTot, filteredTot, residualBBTot] = pwelchFullAndBB(pressure_total, binwidth,fs,P_ref,numofmics, RPM, numblades);

 nooverlap=binwidth*(1-0.5); %50 percent overlap
 nfft=binwidth; %number of points to use in the fft=binwidth.
 dfForPwelch=fs/binwidth; 
 
 [W, ~, ~] = Hann_window(binwidth); %creates hann W sized to binwidth
 
 ef = sum(W.*W)/binwidth; % mean-square of W & 1/Aw = energy_factor. applied automatically within pwelch

%now on the reduced set of data
[pxxRed,fpwelch,PSDRed,PowerSpectrumRed,~,~,oasplRed, spl100Red, spl500Red, spl1000Red, spl10kRed]=pwelchfor1Case(pressure_reduced, W, nooverlap, fs,dfForPwelch,P_ref,numofmics, false);

%% Plots
%plots PSD against fpwelch for unfiltered data

x0=75; %distance from screen when graphing?
y0=75; %same
width=1100; %graph size
height=400; %graph size
graphPos=[x0,y0,width,height];

%all the mics together
createPSDPlotLogScaleSingular(fpwelch, PSDRed, 'Unfiltered PSD Plot (Not downsampled)', 1:numofmics, [0, 30000], [0,70], graphPos, true)
createPSDPlotLogScaleSingular(fpwelch, PSDRed, 'Unfiltered PSD Plot (Not downsampled)', 1:numofmics, [0, 40000], [0,70], graphPos, true)

%all the mics seperately
for i=1:numofmics
    createPSDPlotLogScaleSingular(fpwelch, PSDRed, ['Unfiltered PSD Plot (Not downsampled) for Microphone ',int2str(i)], i, [0, 30000], [0,70], graphPos, true)
    createPSDPlotLogScaleSingular(fpwelch, PSDRed, ['Unfiltered PSD Plot (Not downsampled) for Microphone ',int2str(i)], i, [0, 40000], [0,70], graphPos, true)
end

saveDeletePlots(savePath);

toc
