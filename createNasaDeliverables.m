%Purpose: to create plots NASA asked for (details commented in code)

%Inputs: SPLBPFHarm is the SPLs for the 1-4th harmonics of the BPF. each row a new
%       mic (excluing mic 12), each column is a new harmonic
%       angle is the vector for angle to microphones
%       numofmics is the number of microphones
%       NasaSPL is 1d array where each cell represents SPL in range [500,
%       20000] Hz for each microphone.
%       oaspl is that but for overall sound pressure level, not just after
%       500 Hz.

%Outputs: many graphs

%for testing: createNasaDeliverables(fpwelchTot, pxxTot,SPLForBPFHarmonics, P_ref, angle, numofmics)


% % %function createNasaDeliverables(freqs, pxx, SPLBPFHarm, P_ref, angle, numofmics)
function createNasaDeliverables(SPLBPFHarm, angle, numofmics, NasaSPL, oaspl)

%code from when i had to calculate the [500 to 20000] spl here
% % % freqInterval=[500, 20000]; %interval NASA wants SPL of (in Hz)
% % % startIndex=round(freqInterval(1)/freqs(2))+1; %(fpwelchFiltered(2)=df=difference between frequencies in freqs vector)
% % % endIndex=round(freqInterval(2)/freqs(2))+1;
% % % NasaSPL=zeros(1,(numofmics-1));
% % % oaspl=NasaSPL;
% % % disp(startIndex)
% % % disp(endIndex)
% % % for i=1:(numofmics-1)
% % %     NasaSPL(i)=calculateSPL(freqs,pxx(:,i),P_ref,[startIndex, endIndex]); %SPL in the range NASA desires...
% % %     oaspl(i)=calculateSPL(freqs,pxx(:,i),P_ref,[1, length(freqs)]);
% % % end

%1) SPL in frequency range 500 - 20000 Hz. I think I already did this
%elsewhere, but I'll do it again here just because.
figure()
polarplot(angle(1:(numofmics-1)),NasaSPL(1:(numofmics-1)))
title('500Hz - 20000Hz SPL Directivity')
thetalim([-40 , 30]);
rlim([60 , 80]);


%2) nothing. something already asked for elsewhere...


%3) SPL at BPF, 2nd and 3rd harmonics

% THESE LINES ARE LEFT OVER FROM WHEN I WASN'T CORRECTING EVERYTHING BY
% DISTANCE, SO I DID IT HERE ONLY. NOW, I DO THE CORRECTION IN CurrentCodeforAJHV6 (main)
% % % %need to correct to mic size 1.8956m
% % % SPLBPFHarmCorr=zeros(size(SPLBPFHarm)); %SPL at the harmonics of the BPF, corrected to mic size nasa wants
% % % for i=1:(numofmics-1) %loops through mics except last
% % %     for j=1:length(SPLBPFHarm(1,:)) %loops through each harmonic of the BPF we're looking at 
% % %         SPLBPFHarmCorr(i,j)=SPLBPFHarm(i,j)-abs( 20*log10(micDist(i)/correctDist) );
% % %     end
% % % end

SPLBPFHarmCorr=SPLBPFHarm; %due to the fact that i correct the signal itself, the SPL is the corrected SPL.
figure()
%polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'k' , angle(1:(numofmics-1)),SPL(:,1),'r' , angle(1:(numofmics-1)),SPL(:,2),'g' , angle(1:(numofmics-1)),SPL(:,3),'b' , angle(1:(numofmics-1)),SPL(:,4),'m')
polarplot(angle(1:(numofmics-1)),SPLBPFHarmCorr(:,1),'r' , angle(1:(numofmics-1)),SPLBPFHarmCorr(:,2),'g' , angle(1:(numofmics-1)),SPLBPFHarmCorr(:,3),'b')
thetalim([-40,30]);
rlim([30,80]);
title('BPF Harmonics SPL');
legend('Blade Pass Frequency SPL','2nd Harmonic SPL','3rd Harmonic SPL')

%3.1) same, but including the SPL [500 20000] Hz
figure()
%polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'k' , angle(1:(numofmics-1)),SPL(:,1),'r' , angle(1:(numofmics-1)),SPL(:,2),'g' , angle(1:(numofmics-1)),SPL(:,3),'b' , angle(1:(numofmics-1)),SPL(:,4),'m')
polarplot(angle(1:(numofmics-1)),SPLBPFHarmCorr(:,1),'r' , angle(1:(numofmics-1)),SPLBPFHarmCorr(:,2),'g' , angle(1:(numofmics-1)),SPLBPFHarmCorr(:,3),'b',angle(1:(numofmics-1)),NasaSPL(1:(numofmics-1)),'k')
thetalim([-40,30]);
rlim([30,80]);
title('BPF Harmonics SPL');
legend('Blade Pass Frequency SPL','2nd Harmonic SPL','3rd Harmonic SPL','SPL range [0.5, 20] kHz')

%3.2 and now including the oaspl instead
figure()
%polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'k' , angle(1:(numofmics-1)),SPL(:,1),'r' , angle(1:(numofmics-1)),SPL(:,2),'g' , angle(1:(numofmics-1)),SPL(:,3),'b' , angle(1:(numofmics-1)),SPL(:,4),'m')
polarplot(angle(1:(numofmics-1)),SPLBPFHarmCorr(:,1),'r' , angle(1:(numofmics-1)),SPLBPFHarmCorr(:,2),'g' , angle(1:(numofmics-1)),SPLBPFHarmCorr(:,3),'b',angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'k')
thetalim([-40,30]);
rlim([30,80]);
title('BPF Harmonics SPL');
legend('Blade Pass Frequency SPL','2nd Harmonic SPL','3rd Harmonic SPL','OASPL')

end