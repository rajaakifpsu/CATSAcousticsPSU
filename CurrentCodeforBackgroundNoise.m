%WARNING: user will need to manually update the HPcutoff variable (where
%the highpass filter filters to) when the equivalent variable is changed in
%CurrentCodeforAJHV_. They should match.

%Purpose: to calculate and store acoustics for background noises (multiple
%cases). Then they can be loaded into CurrentCodeforAJHV_ conveniently.

%Input: User will need to input the filepath of whatever they want to
%analyze

%Output: the code will save a mat file containing 3 2d arrays for the PSD of various tests. 
%Columns are mics, rows are frequencies corresponding to the frequencies in fpwelch (a variable used
%in CurrentCodeforAJHV_, frequencies returned by pwelch function). The tests are
%for background with psu, background without psu, and motor noise.

%Misc note: assumes mic sensitivities are the same for all the
%background/motor noise tests
%assuming all the microphones have the same sampling rate



close all ;
clc;

tic


%% Loading in data for performance and acoustics

numofmics = 12;

HPcutoff=100; %frequency value (Hz) to highpass to. 

%background noise with no power supply unit. 
BGnoPSUfilePath='C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_15_21\bg5\acs_data.h5';
BGnoPSU = hdf5info(BGnoPSUfilePath);
%background noise with power supply unit
BGyesPSUfilePath='C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_15_21\bg6\acs_data.h5';
BGyesPSU = hdf5info(BGyesPSUfilePath);
%motor noise
motorNoiseFilePath='C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_15_21\mnbo7\acs_data.h5';
motorNoise = hdf5info(motorNoiseFilePath);



%Reading all the Signals from the h5 file
%to work with the processing techniques
fs = hdf5read(BGnoPSU.GroupHierarchy.Datasets(4)); %sampling rate in Hz. Which is the 4th field
BGnoPSU_acousticsdata_uncorr = hdf5read(BGnoPSU.GroupHierarchy.Datasets(1)); % Acoustics value in voltage.
BGyesPSU_acousticsdata_uncorr = hdf5read(BGyesPSU.GroupHierarchy.Datasets(1)); % Acoustics value in voltage.
motorNoise_acousticsdata_uncorr = hdf5read(motorNoise.GroupHierarchy.Datasets(1)); % Acoustics value in voltage.


%% Downsampling the signals

BGnoPSU_acousticsdata_uncorr = BGnoPSU_acousticsdata_uncorr(1:2:length(BGnoPSU_acousticsdata_uncorr(:,1)),:); %downsampling the signals for computational efficiency
BGyesPSU_acousticsdata_uncorr = BGyesPSU_acousticsdata_uncorr(1:2:length(BGyesPSU_acousticsdata_uncorr(:,1)),:);
motorNoise_acousticsdata_uncorr = motorNoise_acousticsdata_uncorr(1:2:length(motorNoise_acousticsdata_uncorr(:,1)),:);
fs=fs/2; %downsampling the signal for computational efficiency

%% Converting from voltage to Pa with microphone sensitivities

micsensitivities = hdf5read(BGnoPSU.GroupHierarchy.Datasets(5)); % Microphone sensitivities mv/Pa
BGnoPSU_acousticsdata = (1000./(micsensitivities')).*BGnoPSU_acousticsdata_uncorr; % Calibrated Acoustic values%initialize calculations
BGyesPSU_acousticsdata = (1000./(micsensitivities')).*BGyesPSU_acousticsdata_uncorr;
motorNoise_acousticsdata = (1000./(micsensitivities')).*motorNoise_acousticsdata_uncorr;

noofsamples = length(BGnoPSU_acousticsdata); %length of signal

%% Correct the distance from rotor to microphones to be the distance NASA wants.

correctDist=1.8956; %in m (the distance NASA wants)
micDist=[65.19, 62.97, 61.34, 60.34, 60.00, 60.34, 61.34, 62.97, 65.19, 67.93, 71.14, 74.75];%in inches. distance we have to each mic in test.
micDist=convlength(micDist, 'in','m'); %converts micDist to m

for i=1:numofmics
    BGnoPSU_acousticsdata(:,i) = BGnoPSU_acousticsdata(:,i) * (micDist(i)/correctDist); %its just signal multiplied by (dist we have / dist we want)
    BGyesPSU_acousticsdata(:,i) = BGyesPSU_acousticsdata(:,i) * (micDist(i)/correctDist); %its just signal multiplied by (dist we have / dist we want)
    motorNoise_acousticsdata(:,i) = motorNoise_acousticsdata(:,i) * (micDist(i)/correctDist); %its just signal multiplied by (dist we have / dist we want)
end

%% Defining other little variables

N = noofsamples;
dt = 1/fs;    %time increment between each data point
time_array = [0:N-1]*dt;  %get array of times at each point
T = max(time_array); %max time or period

%% take all the pressure time data here

for i = 1:numofmics
    BGnoPSU_pressure_total(:,i)=BGnoPSU_acousticsdata(:,i); 
    BGyesPSU_pressure_total(:,i)=BGyesPSU_acousticsdata(:,i);
    motorNoise_pressure_total(:,i)=motorNoise_acousticsdata(:,i);
end 

%%  taking a section of time in the .h5 file to define pressure_reduced

%turns out this section of code is not necesary since these tests only last
%30s...

%% Correct the DC Offset

for i = 1:numofmics
   
   meanDCoffset(:,i) = mean(BGnoPSU_pressure_total(:,i),1);
   BGnoPSU_pressure_total(:,i) = BGnoPSU_pressure_total(:,i) - meanDCoffset(1,i);
    
   meanDCoffset(:,i) = mean(BGyesPSU_pressure_total(:,i),1);
   BGyesPSU_pressure_total(:,i) = BGyesPSU_pressure_total(:,i) - meanDCoffset(1,i);
   
   meanDCoffset(:,i) = mean(motorNoise_pressure_total(:,i),1);
   motorNoise_pressure_total(:,i) = motorNoise_pressure_total(:,i) - meanDCoffset(1,i);
   
end

%% calculate the OASPL RMS

%THIS MIGHT BE USEFUL, IM NOT SURE. IM NOT USING IT RIGHT NOW, SO ILL
%COMMENT IT OUT FOR NOW, CHANGE IT TO WORK LATER IF NEED BE.

 dbref = 2e-5; % reference pressure in Pa or 2e-5;
 P_ref = dbref;
 G_ref = dbref^2; 
% 
% for k = 1:numofmics
%     rmspress_total(k,:) = rms(pressure_total(1:N,k));
%     OASPL_total_RMS(k,:) = 20*log10(rmspress_total(k,:)/P_ref); %OASPL at each Mic
%     OASPL_reduced_RMS(k,:) = 20*log10(rms(pressure_reduced(1:N_reduced,k))/P_ref); %OASPL reduced at each Mic
%     deltatotalreduced = OASPL_total_RMS - OASPL_reduced_RMS;
% end 


%% Angles of the Microphone Array

angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert angles to radians

%% define RPM
performanceSampRate=20000;
numblades = 2;

%% Highpass filter
BGnoPSU_pressure_total = highpass(BGnoPSU_pressure_total , HPcutoff , fs); 
BGyesPSU_pressure_total = highpass(BGyesPSU_pressure_total , HPcutoff , fs); 
motorNoise_pressure_total = highpass(motorNoise_pressure_total , HPcutoff , fs); 


%% P welch method
 bandwidth = 20; % width of the record in frequency (Hz)
 binwidth = 2^nextpow2(fs/bandwidth); %not necesary for math, but reduces time needed to run code.
 
 nooverlap=binwidth*(1-0.5); %50 percent overlap
 nfft=binwidth; %number of points to use in the fft=binwidth.
 dfForPwelch=fs/binwidth; 
 
 [W, ~, ~] = Hann_window(binwidth); %creates hann W sized to binwidth
 
 %calculated acoustics for reduced background noise (no Power Supply Unit)
 [pxxBGnoPSU,fpwelch,PSDBGnoPSU,PowerSpectrumBGnoPSU,~,~,oasplBGnoPSU, spl500BGnoPSU, spl1000BGnoPSU]=pwelchfor1Case(BGnoPSU_pressure_total, W, nooverlap, fs,dfForPwelch,P_ref,numofmics, false);
 
 %calculated acoustics for reduced background noise (with Power Supply Unit)
 [pxxBGyesPSU,~,PSDBGyesPSU,PowerSpectrumBGyesPSU,~,~,oasplBGyesPSU, spl500BGyesPSU, spl1000BGyesPSU]=pwelchfor1Case(BGyesPSU_pressure_total, W, nooverlap, fs,dfForPwelch,P_ref,numofmics, false);
 
 %calculated acoustics for reduced motor noise
 [pxxMotorNoise,~,PSDMotorNoise,PowerSpectrumMotorNoise,~,~,oasplMotorNoise, spl500MotorNoise, spl1000MotorNoise]=pwelchfor1Case(motorNoise_pressure_total, W, nooverlap, fs,dfForPwelch,P_ref,numofmics, false);
 
 %% Saving Results
 %for now, I'll only save the PSDs... maybe later I'll need more
 
 %they are the reduced versions since only 30s of data
 PSDBGnoPSURed = PSDBGnoPSU;
 PSDBGyesPSURed= PSDBGyesPSU;
 PSDMotorNoiseRed = PSDMotorNoise;
 
 save('NoisePSDs.mat','PSDBGnoPSURed','PSDBGyesPSURed','PSDMotorNoiseRed')
 
