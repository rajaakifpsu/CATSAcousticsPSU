clearvars -except testData
clc
close all
%Version 1-in calculating avErrors, use the mean of the window your in
%calculating avErrors
%Between versions 1 and 2, there is very little difference in the values of
%avErrors (on the order of 10^-4) but v1 is probably more correct? 

%FOR THIS CODE TO WORK, NEED TO DOWNLOAD ALL THE full.txt DOCUMENTS AND PUT
%SAVE THEM TO COMPUTER, NAMED LIKE THE LINE BELLOW THIS ONE.
tests=["tdhcs15.txt","tdhcs1.txt","tdhcs16.txt","tdhcs2.txt","tdhcs3.txt","tdhcs4.txt","tdhcs5.txt","tdhcs6.txt","tdhcs7.txt"];
testnums=[15,1,16,2:7];
sampleRate=20000; %20000 measurements/s

%imports all the data from all the txt files into a 3-dimensional array.
%columns and rows correspond to a single txt file, different pages (3rd index)are
%different txt files.
if ~exist('testData')
    testData=zeros(2380000,27,length(tests));
    for i=1:length(tests)
        testData(:,:,i)=readmatrix(tests(i));
        fprintf('read %f\n',i)
    end
    save('testDataTdhcs','testData','-v7.3');
else
    fprintf('File data already stored\n')
end



%CALCULATING ALL AVERAGE RPMs, related data
%===================================================================
RPMs=zeros(length(tests),7);
%RPMs is matrix of solutions related to RPMs. 
%1st column is names - tdhcs___.
%2nd column is av RPM of motor one. in this case, 0 or NaN since its not
%spinning.
%3rd column in av RPM of motor two.
%Fourth column is desired RPM. (from testMatrixRaja)
%5th column is error to desired RPM
%6th column is measured RPM (From testMatrixRaja)
%7th column is error to measured RPM
%each row is a new test
RPMs(:,1)=testnums';
for i=1:length(tests)
   [RPMs(i,[2,3]),~,~]=getRPMSteady(testData(:,:,i),sampleRate); %gets average RPMs of all tests
end
%these next lines are just copying values from
%TestMatrixRaja_SepEntry_II.xls, comparing to my averages, and displaying
RPMs(:,4)=[1496,1895,1994,2493,2991,3490,3989,4487,4986]';
RPMs(:,5)=RPMs(:,4)-RPMs(:,3);
RPMs(:,6)=[1510,1894,2004,2484,3000,3477,4005,4463,4981,]';
RPMs(:,7)=RPMs(:,6)-RPMs(:,3);
T=array2table(RPMs,'VariableNames',{'Name tdhcs__','m1 RPM', 'm2 RPM', 'desired RPM', 'to desired', 'measured RPM','to measured'})



%CALCULATING AVERAGE FORCES & MOMENTS
%==============================================================
FsMsRPMs=zeros(length(tests),8);
%columns are averages [test# F2x F2y F2z M2x M2y M2z m2RPM] (motor 1 is neglected, not rotating)
%each row is a new test
FsMsRPMs(:,1)=testnums';
for i=1:length(tests)
    temp=getFsMs(testData(:,:,i));
    FsMsRPMs(i,2:7)=temp(7:12); %gets average forces and moments of motor2 for all tests
end
FsMsRPMs(:,8)=RPMs(:,3); %adds in average RPMs of motor 2 
T2=array2table(FsMsRPMs, 'VariableNames',{'Name tdhcs__', 'F2x', 'F2y', 'F2z', 'M2x', 'M2y', 'M2z', 'm2RPM'})
%ave from all 120s


%FINDING ERRORS AND STANDARD DEVIATIONS
%=================================================================
%Strategy: 10s intervals with 5s overlap (10s=200000 measurements)
%For each interval, for all 7 metric (3 forces 3 moments 1 RPM), calculate
%its standard deviation and error to average
%NOTE: I am very unsure of if the way I did the RPM comparison is correct.
%I get an array of the number of measurements between detecting new blades,
%convert that to RPMs, then compare to the average RPM from the whole test.
winlen=10; %length of window (s)
overlap=0.5; %overlap between windows
win=winlen*sampleRate; %length of window (number of measurements)

avErrors=zeros(0);%3D array describing the average error acrued in the 10 window
%each column corresponds error of a different metric (starttime, endtime, 3 forces, 3 moments, 1 RPM)
%each row is a different 10s window
%each page is a different test#.

stanDevs=zeros(0); %same as avErrors, but for standard deviation.

for i=1:length(tests) %loops through all the tests
    startmeasurement=1; %measurment the window starts on
    endmeasurement=win; %measurement the window ends on 
    lm=-1; %incremented variable needed in assigning errorAv to avErrors
    while endmeasurement<length(testData(:,1,1)) %loops through each 10s window.
        lm=lm+1;
        for j=19:24 %loops through all the columns corresponding to forces and moments
            avErrors(lm+1,[1,2],i)=[(startmeasurement-1)/sampleRate,(endmeasurement-1)/sampleRate]; %stores starttime & endtime
            stanDevs(lm+1,[1,2],i)=avErrors(lm+1,[1,2],i);
            errorsum=0; %sum of errors for any variable, eventually stored
            windMean=mean( testData(startmeasurement:endmeasurement,j,i) ); %the mean of the metric in the window
            for k=startmeasurement:endmeasurement %loops through all components of the column in the window
                errorsum=errorsum+abs(testData(k,j,i) - windMean);% adds error of component to errorsum
            end
            errorAv=errorsum/win; %average error for window
            avErrors(startmeasurement-(lm*((win*(1-overlap))-1)),j-16,i)=errorAv;  %stores average error for window
            stanDevs(startmeasurement-(lm*((win*(1-overlap))-1)),j-16,i)=std( testData(startmeasurement:endmeasurement,j,i) ); %calculates and stores standard deviation
        end
        
        [~,AvRPMError,AvRPMStanDev]=getRPMSteady(testData(startmeasurement:endmeasurement,:,i),sampleRate); %passes getRPM function a slice of the matrix corresponding to the desired 10s window
        %note that first index in measurementsBetweenInWindow will not be
        %consistent with the others, since the rotor can start between
        %blades. So start from 2nd index in analysis.
        
        avErrors(startmeasurement-(lm*((win*(1-overlap))-1)),9,i)=AvRPMError(2); %stores motor 2's av RPM error
        stanDevs(startmeasurement-(lm*((win*(1-overlap))-1)),9,i)=AvRPMStanDev(2);
        
        startmeasurement=startmeasurement+(win*(1-overlap));%moves the start and end to next window
        endmeasurement=endmeasurement+(win*(1-overlap));
       
    end
end



%DETERMINING WHAT INTERVALS TO USE
%===============================================================
metricWeights=[0,0,1,0,0,1,1]; %Fx, Fy, Fz, Mx, My, Mz, RPM
%Currently weighted so that Fz, Mz, RPM matter equally in deciding best
%window, nothing else matters.
%[percentAvErrors, idealPercErrorWindows, summedPercentageErrors]=toPercent(avErrors, metricWeights, FsMsRPMs, win);
[idealErrorWindows, summedErrors, errorWindowsSorted]=determineWindow(avErrors, metricWeights);
[idealStanDevWindows, summedStanDevs, stanDevWindowsSorted]=determineWindow(stanDevs, metricWeights);

%IdealWindows=array2table( [testnums',idealPercErrorWindows,idealErrorWindows,idealStanDevWindows],'VariableNames',{'trial number tdhcs__','start (%err)','end (%err)','start (err)','end (err)','start (st dev)','end (st dev)'} )
IdealWindows=array2table( [testnums',idealErrorWindows,idealStanDevWindows],'VariableNames',{'trial number tdhcs__','start (err)','end (err)','start (st dev)','end (st dev)'} )

DIVINEwindows=[testnums',zeros(length(tests),2)];

for i=1:length(tests) %loops through the tests
    flag=true;
    for j=1:length(errorWindowsSorted(:,1,1)) %loops through the windows. Note since it'll be working with the sorted arrays, its going from best to worst window, not chronological
        %compares starttimes, to see if they are in the top j in both
        %sorted lists (errorWindowsSorted and StanDevWindowsSorted)
        temp=intersect(errorWindowsSorted(1:j,1,i) , stanDevWindowsSorted(1:j,1,i));
        if (flag && (length(temp)==1) ) %if the one of jth best windows in the two lists are the same, thats the best window
            DIVINEwindows(i,[2,3])=[temp, temp+winlen];
            flag=false;
        elseif (flag && (length(temp)==2)) %if two intersections are uncovered, need to choose one. 
            temp1Score=find(errorWindowsSorted(:,1,i)==temp(1)) + find(stanDevWindowsSorted(:,1,i)==temp(1)); %score is sum of indices in sorted list
            temp2Score=find(errorWindowsSorted(:,1,i)==temp(2)) + find(stanDevWindowsSorted(:,1,i)==temp(2)); %I just now realize I could have done this whole part this way... no turning back now!
            DIVINEwindows(i,[2,3])=[temp(1),temp(1)+winlen];
            if temp2Score>temp1Score
                DIVINEwindows(i,[2,3])=[temp(2),temp(2)+winlen];
            end
            flag=false;
        end
    end
end
save('DIVINEwindowsTdhcs.mat','DIVINEwindows')
TDivineWindows=array2table(DIVINEwindows, 'VariableNames',{'test name tdhcs__','start time','end time'})





%Plotting stuff
%====================================================
% temp=[13.4,13.4,13.4,13.4,13.4,13.4,13.4,13.4,13.4]; %tempertures (C) for the tdhcs tests, in order [15,16,2:7]. Manually copied from TestMatrixRaja_SepEntry_II
% press=[976,976,976,976,976,976,976,976,976];%same but for pressures (mbar)
% createFnMProfiles(testData, DIVINEwindows, win, sampleRate, 'TDHCS ', 0, temp,press,false,false)
% createFnMProfiles(testData, DIVINEwindows, win, sampleRate, 'TDHCS ', 0, temp,press,true,false) %the way createFnMProfiles is set up, the figures generated by these two will be superimposed.






%CONCLUSION: This method works pretty good for most tests. 
%Test 3 it doesn't work on too too well, but I don't think there's a better
%10s window to choose from. No overlap between ranking of error and st dev
%for the top 4-5 intervals.
%I also did a little work with percent error, but the intervals it said
%were good were not the intervals that they other 2 methods said were good,
%so i abandoned it.

%TO DO WITH GOOD CODE: insert the RPMs from getRPM (measured RPM)(for the 10s window decided on here.)
%take data from DIVINEwindows and put it into profile for all forces and moments. Have 2 versions of both: corrected for density and uncorrected
%for density 
%(instead of the desired RPM, use the actual measured RPM for the
%window)(just pass getRPM the 10s window of data I want)
%also include the SEM-standard error measurment
%Then we can move on to processing the case with 2 blades.















%Junk code, disregard. Saving here just in case.
% clearvars -except testData
% clc
% %Version 1-in calculating avErrors, use the mean of the window your in
% %calculating avErrors
% %Between versions 1 and 2, there is very little difference in the values of
% %avErrors (on the order of 10^-4) but v1 is probably more correct? 
% 
% %FOR THIS CODE TO WORK, NEED TO DOWNLOAD ALL THE full.txt DOCUMENTS AND PUT
% %SAVE THEM TO COMPUTER, NAMED LIKE THE LINE BELLOW THIS ONE.
% tests=["tdhcs15.txt","tdhcs1.txt","tdhcs16.txt","tdhcs2.txt","tdhcs3.txt","tdhcs4.txt","tdhcs5.txt","tdhcs6.txt","tdhcs7.txt"];
% testnums=[15,1,16,2:7];
% sampleRate=20000; %20000 measurements/s
% 
% %imports all the data from all the txt files into a 3-dimensional array.
% %columns and rows correspond to a single txt file, different pages (3rd index)are
% %different txt files.
% if ~exist('testData')
%     testData=zeros(2380000,27,length(tests));
%     for i=1:length(tests)
%         testData(:,:,i)=readmatrix(tests(i));
%         fprintf('read %f\n',i)
%     end
% else
%     fprintf('File data already stored\n')
% end
% 
% 
% 
% %CALCULATING ALL AVERAGE RPMs, related data
% %===================================================================
% RPMs=zeros(length(tests),7);
% %RPMs is matrix of solutions related to RPMs. 
% %1st column is names - tdhcs___.
% %2nd column is av RPM of motor one. in this case, 0 or NaN since its not
% %spinning.
% %3rd column in av RPM of motor two.
% %Fourth column is desired RPM. (from testMatrixRaja)
% %5th column is error to desired RPM
% %6th column is measured RPM (From testMatrixRaja)
% %7th column is error to measured RPM
% %each row is a new test
% RPMs(:,1)=testnums';
% for i=1:length(tests)
%    [RPMs(i,[2,3]),~,~]=getRPMSteady(testData(:,:,i),sampleRate); %gets average RPMs of all tests
% end
% %these next lines are just copying values from
% %TestMatrixRaja_SepEntry_II.xls, comparing to my averages, and displaying
% RPMs(:,4)=[1895,2493,2991,3490,3989,4487,4986,1496,1994]';
% RPMs(:,5)=RPMs(:,4)-RPMs(:,3);
% RPMs(:,6)=[1894,2484,3000,3477,4005,4463,4981,1510,2004]';
% RPMs(:,7)=RPMs(:,6)-RPMs(:,3);
% T=array2table(RPMs,'VariableNames',{'Name tdhcs__','m1 RPM', 'm2 RPM', 'desired RPM', 'to desired', 'measured RPM','to measured'})
% 
% 
% 
% %CALCULATING AVERAGE FORCES & MOMENTS
% %make this a function accepting 1 or 2 for motor (eventually for 2-motor case)
% %==============================================================
% FsMsRPMs=zeros(length(tests),8);
% %columns are averages [test# F2x F2y F2z M2x M2y M2z m2RPM] (motor 1 is neglected, not rotating)
% %each row is a new test
% FsMsRPMs(:,1)=testnums';
% for i=1:length(tests)
%    temp=getFsMs(testData(:,:,i));
%    FsMsRPMs(i,2:7)=temp(7:12); %gets average forces and moments of motor2 for all tests
% end
% FsMsRPMs(:,8)=RPMs(:,3); %adds in average RPMs of motor 2 
% T2=array2table(FsMsRPMs, 'VariableNames',{'Name tdhcs__', 'F2x', 'F2y', 'F2z', 'M2x', 'M2y', 'M2z', 'm2RPM'})
% %ave from all 120s
% 
% 
% %FINDING ERRORS AND STANDARD DEVIATIONS
% %=================================================================
% %Strategy: 10s intervals with 5s overlap (10s=200000 measurements)
% %For each interval, for all 7 metric (3 forces 3 moments 1 RPM), calculate
% %its standard deviation and error to average
% %NOTE: I am very unsure of if the way I did the RPM comparison is correct.
% %I get an array of the number of measurements between detecting new blades,
% %convert that to RPMs, then compare to the average RPM from the whole test.
% winlen=10; %length of window (s)
% overlap=0.5; %overlap between windows
% win=winlen*sampleRate; %length of window (number of measurements)
% 
% avErrors=zeros(0);%3D array describing the average error acrued in the 10 window
% %each column corresponds error of a different metric (starttime, endtime, 3 forces in motor1, 3 moments motor1, 3 forces motor2, 3 moments motor2, 1 RPM motor1, 1 RPM motor2)
% %each row is a different 10s window
% %each page is a different test#.
% 
% stanDevs=zeros(0); %same as avErrors, but for standard deviation.
% 
% for i=1:length(tests) %loops through all the tests
%     startmeasurement=1; %measurment the window starts on
%     endmeasurement=win; %measurement the window ends on 
%     lm=-1; %incremented variable needed in assigning errorAv to avErrors
%     while endmeasurement<length(testData(:,1,1)) %loops through each 10s window.
%         lm=lm+1;
%         for j=19:24 %loops through all the columns corresponding to forces and moments
%             avErrors(lm+1,[1,2],i)=[(startmeasurement-1)/sampleRate,(endmeasurement-1)/sampleRate]; %stores starttime & endtime
%             stanDevs(lm+1,[1,2],i)=avErrors(lm+1,[1,2],i);
%             errorsum=0; %sum of errors for any variable, eventually stored
%             windMean=mean( testData(startmeasurement:endmeasurement,j,i) ); %the mean of the metric in the window
%             for k=startmeasurement:endmeasurement %loops through all components of the column in the window
%                 errorsum=errorsum+abs(testData(k,j,i) - windMean);% adds error of component to errorsum
%             end
%             errorAv=errorsum/win; %average error for window
%             avErrors(startmeasurement-(lm*((win*(1-overlap))-1)),j-16,i)=errorAv;  %stores average error for window
%             stanDevs(startmeasurement-(lm*((win*(1-overlap))-1)),j-16,i)=std( testData(startmeasurement:endmeasurement,j,i) ); %calculates and stores standard deviation
%         end
%         
%         [~,~, measurementsBetweenInWindow]=getRPMSteady(testData(startmeasurement:endmeasurement,:,i),sampleRate); %passes getRPM function a slice of the matrix corresponding to the desired 10s window
%         %note that first index in measurementsBetweenInWindow will not be
%         %consistent with the others, since the rotor can start between
%         %blades. So start from 2nd index in analysis.
%         RPMsInWindow=(measurementsBetweenInWindow(1,:)/(sampleRate*60)).^-1; %converts to RPMs for each revolution in 10s window
%         RPMErrorSum=0; % sum of error in RPMs
%         tempmean=mean(RPMsInWindow(2:length(RPMsInWindow)));
%         for n=2:length(measurementsBetweenInWindow)
%             RPMErrorSum=RPMErrorSum+abs(RPMsInWindow(n)-tempmean);%sums up error between RPMs in window and average RPM of window
%         end
%         AvRPMError=RPMErrorSum/(length(RPMsInWindow)-1); %Average error in RPMs
%         avErrors(startmeasurement-(lm*((win*(1-overlap))-1)),15,i)=AvRPMError; %stores average error in RPMs
%         stanDevs(startmeasurement-(lm*((win*(1-overlap))-1)),15,i)=std( RPMsInWindow(2:length(RPMsInWindow)) ); %calculates and stores standard deviation of RPMs
% 
%     end
% end
% 
% 
% 
% %DETERMINING WHAT INTERVALS TO USE
% %===============================================================
% metricWeights=[0,0,1,0,0,1,1]; %Fx, Fy, Fz, Mx, My, Mz, RPM
% %Currently weighted so that Fz, Mz, RPM matter equally in deciding best
% %window, nothing else matters.
% %[percentAvErrors, idealPercErrorWindows, summedPercentageErrors]=toPercent(avErrors, metricWeights, FsMsRPMs, win);
% [idealErrorWindows, summedErrors, errorWindowsSorted]=determineWindow(avErrors, metricWeights);
% [idealStanDevWindows, summedStanDevs, stanDevWindowsSorted]=determineWindow(stanDevs, metricWeights);
% 
% %IdealWindows=array2table( [testnums',idealPercErrorWindows,idealErrorWindows,idealStanDevWindows],'VariableNames',{'trial number tdhcs__','start (%err)','end (%err)','start (err)','end (err)','start (st dev)','end (st dev)'} )
% IdealWindows=array2table( [testnums',idealErrorWindows,idealStanDevWindows],'VariableNames',{'trial number tdhcs__','start (err)','end (err)','start (st dev)','end (st dev)'} )
% 
% DIVINEwindows=[testnums',zeros(length(tests),2)];
% 
% for i=1:length(tests) %loops through the tests
%     flag=true;
%     for j=1:length(errorWindowsSorted(:,1,1)) %loops through the windows. Note since it'll be working with the sorted arrays, its going from best to worst window, not chronological
%         %compares starttimes, to see if they are in the top j in both
%         %sorted lists (errorWindowsSorted and StanDevWindowsSorted)
%         temp=intersect(errorWindowsSorted(1:j,1,i) , stanDevWindowsSorted(1:j,1,i));
%         if (flag && (length(temp)==1) ) %if the one of jth best windows in the two lists are the same, thats the best window
%             DIVINEwindows(i,[2,3])=[temp, temp+winlen];
%             flag=false;
%         elseif (flag && (length(temp)==2)) %if two intersections are uncovered, need to choose one. 
%             temp1Score=find(errorWindowsSorted(:,1,i)==temp(1)) + find(stanDevWindowsSorted(:,1,i)==temp(1)); %score is sum of indices in sorted list
%             temp2Score=find(errorWindowsSorted(:,1,i)==temp(2)) + find(stanDevWindowsSorted(:,1,i)==temp(2)); %I just now realize I could have done this whole part this way... no turning back now!
%             DIVINEwindows(i,[2,3])=[temp(1),temp(1)+winlen];
%             if temp2Score>temp1Score
%                 DIVINEwindows(i,[2,3])=[temp(2),temp(2)+winlen];
%             end
%             flag=false;
%         end
%     end
% end
% save('DIVINEwindows.mat','DIVINEwindows')
% TDivineWindows=array2table(DIVINEwindows, 'VariableNames',{'test name tdhcs__','start time','end time'})
% 
% %CONCLUSION: This method works pretty good for most tests. 
% %Test 3 it doesn't work on too too well, but I don't think there's a better
% %10s window to choose from. No overlap between ranking of error and st dev
% %for the top 4-5 intervals.
% %I also did a little work with percent error, but the intervals it said
% %were good were not the intervals that they other 2 methods said were good,
% %so i abandoned it.
% 
% %TO DO WITH GOOD CODE: insert the RPMs from getRPM (measured RPM)(for the 10s window decided on here.)
% %take data from DIVINEwindows and put it into profile for all forces and moments. Have 2 versions of both: corrected for density and uncorrected
% %for density 
% %(instead of the desired RPM, use the actual measured RPM for the
% %window)(just pass getRPM the 10s window of data I want)
% %also include the SEM-standard error measurment
% %Then we can move on to processing the case with 2 blades.