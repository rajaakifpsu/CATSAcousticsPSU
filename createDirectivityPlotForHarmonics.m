%Purpose: to create the directivity plots for the harmonics. It also
%       returns the SPL of the BPF harmonics.

%Inputs: pxx is matrix representing pxx for all mics. each row is a freq,
%       each column is a microphone.
%       oaspl is the vector representing the overall sound pressure level of each mic.
%       fpwelch is the vector of frequencies corresponding to the above
%       variables.
%       RPMmotor is only the RPM of motor2, the one that is spinning, NOT [motor1 motor2]
%       numblades is the number of blades
%       numofmics is the number of microphones
%       P_ref is the reference pressure
%       graphtitle is title of the graph

%Outputs: SPLforHarm is a 2d array representing the Sound pressure level of
%       different harmonics. Each row is a mic (excluding last mic since it is a 
%       different size.) Each column is a harmonic, from 1 to 4.
%       A directivity plot for OASPL, BPF, harmonic 2, harmonic 3. Note
%       that MATLAB doesnt allow for axis titles, so they will need to be
%       added in manually.
%       A bunch of plots of pxx with the frequencies for each BPF harmonic
%       highlighted in other colors.

%for testing:
%createDirectivityPlotForHarmonics(pxxTot,oasplTot,fpwelchTot,RPM(2),numblades,numofmics,P_ref,'BPF Harmonics SPL Directivity Plot')

function SPLforHarm=createDirectivityPlotForHarmonics(pxx,oaspl,fpwelch,RPMmotor,numblades,numofmics,P_ref,graphtitle)

dfpwelch=fpwelch(2)-fpwelch(1); 

angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert into radians (only for plotting)

BPFs=getBPFHarmonics(numblades,RPMmotor);
%first column is 1st harmonic, 2nd is 2nd harmonic, 3rd is 3rd harmonic,
%4th is 4th harmonics

SPLforHarm=zeros(numofmics-1,4); %SPL for 1st 2nd 3rd 4th harmonic. Each row different mic.

freqIndexLims=getBPFfreqIndexLims(BPFs, dfpwelch, fpwelch, pxx(:,1)); %the BPF and harmonic tones will occur at same indices for all microphones.

for h=1:(numofmics-1) %loops through mics, ignoring last mic since its different size.
    %Doing reiman sum to aproximate the SPL of the tones corresponding to
    %frequencies of BPF and its harmonics
    for i=1:length(BPFs) %loops through the harmonics
        %shoudl be freqIndexLims(2,i)  bellow ??? no i dont think so...
        SPLforHarm(h,i)=calculateSPL(fpwelch,pxx(:,h),P_ref,[freqIndexLims(1,i) , freqIndexLims(2,i)]); %calculates SPL for each harmonic of each mic once looped through a bunch of times
        %old code: still works but in function now.
% %         rsum=0;
% %         for j=freqIndexLims(1,i):(freqIndexLims(2,i)-1) %loops through frequencies corresponding to the harmonic
% %             rsum=rsum+(.5 * (fpwelch(j+1)-fpwelch(j)) * (pxx(j,h)+pxx(j+1,h)) );
% %         end
% %         SPL(h,i)=10*log10(rsum/(P_ref^2)); %not sure if i should be using this since it matches whats used in oaspl...
% %         %SPL(h,i)=rsum;
    end
    
end

figure()
%polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'k' , angle(1:(numofmics-1)),SPL(:,1),'r' , angle(1:(numofmics-1)),SPL(:,2),'g' , angle(1:(numofmics-1)),SPL(:,3),'b' , angle(1:(numofmics-1)),SPL(:,4),'m')
polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'k' , angle(1:(numofmics-1)),SPLforHarm(:,1),'r' , angle(1:(numofmics-1)),SPLforHarm(:,2),'g' , angle(1:(numofmics-1)),SPLforHarm(:,3),'b')
thetalim([-40,30]);
rlim([30,80]);
title(graphtitle);
legend('OASPL','Blade Pass Frequency SPL','2nd Harmonic SPL','3rd Harmonic SPL')

%%%disp(pxx(1:10,:))

%plot useful for seeing what is being counted as the frequency ranges of
%the different BPF harmonics.
for i=1:numofmics
    figure()
    hold on;
    plot(fpwelch, pxx(:,i) , 'k')
    plot(fpwelch(freqIndexLims(1,1):freqIndexLims(2,1)), pxx(freqIndexLims(1,1):freqIndexLims(2,1),i), 'r')
    plot(fpwelch(freqIndexLims(1,2):freqIndexLims(2,2)), pxx(freqIndexLims(1,2):freqIndexLims(2,2),i), 'g')
    plot(fpwelch(freqIndexLims(1,3):freqIndexLims(2,3)), pxx(freqIndexLims(1,3):freqIndexLims(2,3),i), 'b')
    plot(fpwelch(freqIndexLims(1,4):freqIndexLims(2,4)), pxx(freqIndexLims(1,4):freqIndexLims(2,4),i), 'm')
    xlim([0,5000]);
    ylim([0,10^(-4)])
    ylabel('Pxx (Pa^2/Hz)')
    xlabel('frequency (Hz)')
    title(['Pxx v frequency (mic ' , int2str(i),')'])
    hold off
end

end