%purpose: to create force and moment profiles (graphing those vs RPM).
%input :tData=all the data from all the Full.txt's
%       windows=the best windows to base everything else off. DIVINEwindows
%       in main code.
%       numMeas=number of measurements per window
%       sampRate= sample rate
%       testTitle=used in the titles of graphs
%       figOffset=to determine whixh figure it starts at.
%       temperture=array of tempertures from each test to graphed, manually copied
%       from TestMatrixRaja_SepEntry_II
%       pressures=same but for pressure
%       correction= boolean for whether or not to do density correction
%       isTwoMotors= true/false for if there's two motors running,
%       determines which get figures get saved.
         

function createFnMProfiles(tData, windows, numMeas, sampRate, testTitle, figOffset, tempertures,pressures,correction, isTwoMotors)

%FINDING DATA I WANT BASED OFF WINDOWS
%===================================================
goodData=zeros(numMeas,length(tData(1,:,1)),length(tData(1,1,:)));
for i=1:length(windows(:,1)) %loops through all tests
    startMeas=1+(windows(i,2)*sampRate);
    endMeas=windows(i,3)*sampRate;
    goodData(:,:,i)=tData(startMeas:endMeas,:,i);
end



%CORRECTING FOR DENSITY
%===================================================
linetype=':o';
if (correction)
    linetype='-o';
    R = 287.058; %J/(kg*K);
    Ts = tempertures+273.15; %converting temperture C-->K
    Ps = pressures*100; %converting mbar-->Pa
    density_psu = Ps./(R.*Ts);  %lower than the std day. at 1148ft or 350m. 
    density_stdday = 1.225; %kg/m3
    for i=1:length(goodData(1,1,:)) %loops through each test
        goodData(:,13:24,i)=goodData(:,13:24,i).*density_stdday./density_psu(i); %does density correction on forces and moments
    end
end



%Getting averages and SEMs
%==================================================
FMRPMWindowed=zeros(length(tData(1,1,:)),14); %each column is average of a different metric 3Fs, 3Ms, 3Fs, 3Ms, 2RPMs. Each row is a test
stanDevs=zeros(length(tData(1,1,:)),12); %same but standard deviation.
SEM=stanDevs; %standard error margin
for i=1:length(windows(:,1)) %loops through all tests
    FMRPMWindowed(i,1:12)=getFsMs(goodData(:,:,i));
    [RPMs,~,~]=getRPMSteady(goodData(:,:,i),sampRate);
    FMRPMWindowed(i,[13,14])=RPMs;
    for j=1:12 %loops through all forces and moments
        stanDevs(i,j)=std(goodData(:,j,i));
        SEM(i,j)=stanDevs(i,j)/(sqrt(numMeas));
    end
end

totalGoodData=zeros(size(goodData));
totalGoodData(:,13:18,:)=(goodData(:,13:18,:)+goodData(:,19:24,:)); %adds motor 1 and 2. 0s elsewhere so it works with other functions. last 2 columns are RPMs
totalFMRPMWindowed=zeros(length(tData(1,1,:)),7); %each column is average of a different metric 3Fs, 3Ms, RPM. Each row is a test
totalStanDevs=zeros(length(tData(1,1,:)),12); %same but standard deviation.
totalSEM=totalStanDevs; %standard error margin
for i=1:length(windows(:,1)) %loops through all tests
    temporary=getFsMs(totalGoodData(:,:,i));
    totalFMRPMWindowed(i,1:6)=temporary(1:6);
    totalFMRPMWindowed(i,7)=(FMRPMWindowed(i,13)+FMRPMWindowed(i,14))/2; %total RPM is average of both motors??
    for j=1:6 %loops through all forces and moments
        totalStanDevs(i,j)=std(totalGoodData(:,j,i));
        totalSEM(i,j)=totalStanDevs(i,j)/(sqrt(numMeas));
    end
end

%Plotting the stuff
%===================================================
yaxis=["Fx","Fy","Fz","Mx","My","Mz"];
titles=["Fx v RPM" , "Fy v RPM","Fz v RPM","Mx v RPM","My v RPM","Mz v RPM"];
motornum=["Motor 1 ", "Motor 2 ", "Total "];
for i=1:6 %loops through the 6 forces/moments
    temptitle=append(testTitle,' ',motornum(1),titles(i));
    figure(figOffset+(i*3)-2)
    grid on;
    hold on
    errorbar( FMRPMWindowed(:,13) , FMRPMWindowed(:,i) , SEM(:,i),strcat('b',linetype)) %plots the f/m vs RPM for motor1
    xlabel('RPM')
    ylabel(yaxis(i))
    title(temptitle)
    if (correction && isTwoMotors)
        saveas(gcf , strcat('Profiles/',temptitle,'.pdf'))
    end
    
    temptitle=append(testTitle,' ',motornum(2),titles(i));
    figure(figOffset+(i*3)-1)
    grid on;
    hold on;
    errorbar( FMRPMWindowed(:,14) , FMRPMWindowed(:,i+6), SEM(:,i+6), strcat('r',linetype)) %plots f/m vs RPM for motor2
    xlabel('RPM')
    ylabel(yaxis(i))
    title(temptitle)
    if (correction)
        saveas(gcf , strcat('Profiles/',temptitle,'.pdf'))
    end
    
    temptitle=append(testTitle,' ',motornum(3),titles(i));
    figure(figOffset+(i*3))
    grid on;
    hold on;
    errorbar( totalFMRPMWindowed(:,7) , totalFMRPMWindowed(:,i) , totalSEM(:,i), strcat('g',linetype)) %plots f/m vs RPM for total
    xlabel('RPM')
    ylabel(yaxis(i))
    title(temptitle)
    if (correction && isTwoMotors)
        saveas(gcf , strcat('Profiles/',temptitle,'.pdf'))
    end

%     figure(i)
%     grid on;
%     hold on;
%     errorbar( FMRPMWindowed(:,13) , FMRPMWindowed(:,i) , SEM(:,i),strcat('b',linetype)) %plots the f/m vs RPM for motor1
%     errorbar( FMRPMWindowed(:,14) , FMRPMWindowed(:,i+6), SEM(:,i+6), strcat('r',linetype)) %plots f/m vs RPM for motor2
%     errorbar( totalFMRPMWindowed(:,7) , totalFMRPMWindowed(:,i) , totalSEM(:,i), strcat('g',linetype)) %plots f/m vs RPM for total
%     xlabel('RPM')
%     ylabel(yaxis(i))
%     legend('Motor 1','Motor 2','Total')
%     title(strcat(testTitle,' ',titles(i)))
%     hold off;
end

%%%% createFnMProfiles(testData(:,:,1:8),DIVINEwindows(1:8,:),win,sampleRate,'bahhhhhh',[2 4 10],[1,1,1],false) %%%used for testing purposes on testData generated by
%AnalysingCCWSteady2plus2V1.m
end