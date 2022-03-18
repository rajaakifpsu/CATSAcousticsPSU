%inputs:
%       RPMs1 from getRPMRamp
%       RPMs2 from getRPMRamp
%       fig can be 1xwhatever or 2xwhatever depending on if its single-or
%       double motor, (whatever=number of tests)
%       chunkPortion, overlap,used in smoothInTimeDomain
%       titles used fot titles
%       groupLen = number of tests in a group
%xlimit here is only 2D, a single page of xlimit in createGraph9

function createGraph4(RPMs1,RPMs2,tData,tDataCorr, xlimit, sampRate, groupLen, chunkPortion, overlap, titles,fig)

%for motor 2 (the one that always runs)
%for i=1:1 %only does first test. for speed while testing code
for i=1:length(tData(1,1,:)) %loops through all the tests
    tData(:,1,i)=tData(:,1,i)-tData(1,1,i); %makes it so the test starts at time=0
    tDataCorr(:,1,i)=tData(:,1,i);
    stop=300000;
    toDisplay=[1:stop]; %only need to plot part of data, otherwise computer dies.
    fprintf('%f \n',i)
    
    figure(fig(1,i))
    
    subplot(3,1,1)
    hold on;
    grid on;
    %plot(tData(toDisplay,1,i) , tData(toDisplay,21,i),'-r') %plots Fz vs time
    plot(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,21,i),'-r') %Fz vs time (density-corrected)
    xlabel('Time (s)');
    ylabel('Thrust (N)');
    xlim(xlimit(i,:));  
    FzFixed=smoothInTimeDomain(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,21,i),chunkPortion, overlap); %creates a data for a smoothed version of Fz plot
    plot(FzFixed(2,:),FzFixed(1,:),'-k') %plots it
    hold off;
    
    subplot(3,1,2)
    hold on;
    grid on;
    %plot(tData(toDisplay,1,i) , tData(toDisplay,24,i),'-r') %plots Mz vs time
    plot(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,24,i),'-r') %Mz v time (density-corrected)
    xlabel('Time (s)');
    ylabel('Torque (Nm)');
    xlim(xlimit(i,:));
    MzFixed=smoothInTimeDomain(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,24,i) ,chunkPortion, overlap); %creates data for smoothed version of Mz 
    plot(MzFixed(2,:) , MzFixed(1,:),'-k')
    hold off;
    
    subplot(3,1,3)
    hold on;
    grid on;
    stop=find( RPMs2(1,:,i)==0 , 1,'first')-1; %looks for when RPM goes to zero in the matrix and stops graphing after that point (since thats when the data stopped being collected. after that is just a bunch of 0s that need to be there since each test will have a different number of RPMs recorded)
    if isempty(stop)
        stop=length(RPMs2(1,:,i)); %if the RPM never goes to zero (has recordings for all positions in the matrix), then stop at the end
    end
    plot(RPMs2(2,1:stop,i) , RPMs2(1,1:stop,i),'-k'); %plots RPMs against time
    xlabel('Time (s)')
    ylabel('RPM')
    xlim(xlimit(i,:))
    ylim([0,5050]);
    hold off;

    set(gcf, 'position',[10,10,800,1200] ) 
    sgtitle(titles(i,2));
    
end



if ~isempty(RPMs1) %if motor 1 is aslo rotating
    
    for i=1:length(tData(1,1,:)) %loops through all the tests
        tData(:,1,i)=tData(:,1,i)-tData(1,1,i); %makes it so the test starts at time=0
        tDataCorr(:,1,i)=tData(:,1,i);
        stop=300000;
        toDisplay=[1:stop]; %only need to plot part of data, otherwise computer dies.
        fprintf('%f \n',i)

        figure(fig(2,i))

        subplot(3,1,1)
        hold on;
        grid on;
        %plot(tData(toDisplay,1,i) , tData(toDisplay,21,i),'-r') %plots Fz vs time
        plot(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,15,i),'-r') %Fz vs time (density-corrected)
        xlabel('Time (s)');
        ylabel('Thrust (N)');
        xlim(xlimit(i+groupLen,:));  
        FzFixed=smoothInTimeDomain(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,15,i),chunkPortion, overlap); %creates a data for a smoothed version of Fz plot
        plot(FzFixed(2,:),FzFixed(1,:),'-k') %plots it
        hold off;

        subplot(3,1,2)
        hold on;
        grid on;
        %plot(tData(toDisplay,1,i) , tData(toDisplay,24,i),'-r') %plots Mz vs time
        plot(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,18,i),'-r') %Mz v time (density-corrected)
        xlabel('Time (s)');
        ylabel('Torque (Nm)');
        xlim(xlimit(i+groupLen,:));
        MzFixed=smoothInTimeDomain(tDataCorr(toDisplay,1,i) , tDataCorr(toDisplay,18,i),chunkPortion, overlap); %creates data for smoothed version of Mz 
        plot(MzFixed(2,:) , MzFixed(1,:),'-k')
        hold off;

        subplot(3,1,3)
        hold on;
        grid on;
        stop=find( RPMs1(1,:,i)==0 , 1,'first')-1; %looks for when RPM goes to zero in the matrix and stops graphing after that point (since thats when the data stopped being collected. after that is just a bunch of 0s that need to be there since each test will have a different number of RPMs recorded)
        if isempty(stop)
            stop=length(RPMs1(1,:,i)); %if the RPM never goes to zero (has recordings for all positions in the matrix), then stop at the end
        end
        plot(RPMs1(2,1:stop,i) , RPMs1(1,1:stop,i),'-k'); %plots RPMs against time
        xlabel('Time (s)')
        ylabel('RPM')
        xlim(xlimit(i+groupLen,:)) 
        ylim([0,5050]);
        hold off;

        set(gcf, 'position',[10,10,800,1200] )
        sgtitle(titles(i,1))
    end
    
    
end



end





%5th attempt: now I'll just try averaging in time domain. I'll take the
    %part of data i want to analyze, break it up into a bunch of Hann
    %windows w 50% overlap. For each one, I'll find weighted average of Fz
    %(weighted by Hann coeficients), then plot them at time in center of Hann
    %Window.
% % %     FzToAnalyze=tDataCorr(toDisplay,21,i);
% % %     MeasLength=length(FzToAnalyze); %length of data to analyze (in measurments)
% % %     chunkSize=MeasLength/1000;
% % %     overlap=0.5; %number of overlapping measurments
% % %     increment=chunkSize*(1-overlap);
% % %     loop=1;
% % %     startIndex=1;
% % %     endIndex=chunkSize;
% % %     FzFixed=zeros(0); %first row is Fz, 2nd row is time
% % %     while endIndex<MeasLength %while we havent surpassed the end of data we want analyzed, loops through making new windows
% % %         thisHann=hann(chunkSize);
% % %         thisWindow=FzToAnalyze(startIndex:endIndex).*thisHann; %applies hann window to Fz
% % %         FzAv=sum(thisWindow)/sum(thisHann); %weighted average of Fz in the window
% % %         FzFixed(1,loop)=FzAv; %records the Fz averaged (weighted with hann window)
% % %         
% % %         startTime=tDataCorr(startIndex,1,i);
% % %         endTime=tDataCorr(endIndex,1,i);
% % %         FzFixed(2,loop)=(startTime+endTime)/2; %records time in middle of hann window
% % %         
% % %         startIndex=startIndex+increment;
% % %         endIndex=endIndex+increment;
% % %         loop=loop+1;
% % %     end




 %okay, now trying to determine where it the ramp starts and stops.
%     startRampTime=0;
%     endRampTime=0;
%     startRampIndexRPM=0;
%     cont=true;
%     for j=1:(length(RPMs2(1,:,i))-1) %loops through all RPM measurments
%         changeinRPM=RPMs2(1,j+1,i)-RPMs2(1,j,i);
%         if changeinRPM>=0 && cont %if the change in RPM is positive, then ramp has started
%             startRampIndexRPM=j;
%             fprintf('startramp index %i\n',startRampIndexRPM)
%             startRampTime=RPMs2(2,j,i);
%             cont=false;
%         end
%     end
%     cont=true;
%     for j=startRampIndexRPM:(length(RPMs2(1,:,i))-4) %loops through remaining RPM measurments
%         changeinRPM=RPMs2(1,(j+1):(j+4),i)-RPMs2(1,j,i);
%         if (~isempty(find(changeinRPM<0)) && cont) %if the change between this RPM and any of the next 4 is negative, then ramp ends
%             endRampTime=RPMs2(2,j,i);
%             fprintf('endRampTime: %f\n',endRampTime)
%             cont=false;
%         end
%     end
%     startRampIndex=1+floor(startRampTime*sampRate); %index where ramp starts (in Fz, Mz measurments). floor needed since startRampTime will not give an integer since its inbetween measuremnets
%     endRampIndex=1+floor(endRampTime*sampRate);
%     fprintf('%f Ramp from %f to %f',i,startRampTime,endRampTime)
%     %we don't care about what happens before ramp since rotors not even
%     %supposed to be spinning then.
%     chunkSize=1000; %1/20th a second
%     numOverlap=0.5*chunkSize; %number of overlapping measurments
% % % %     startIndex=startRampIndex;
% % % %     endIndex=startIndex+chunkSize;
% % % %     loop=1;
% % % %     Y=zeros(0);%ffts. each row is a chunk, each column a frequency.
% % % %     while startIndex<endRampIndex %loops through all chunks. Note that startIndex is incremented until it exceeds endRampIndex.
% % % %         %to do: hann window
% % % %         Fz=tDataCorr(startIndex:endIndex , 21, i); %the Fz in the chunk
% % % %         L=length(Fz);%length of signal (in measurements)
% % % %         YinChunk=fft(Fz,L); %fourier transform of the Fz
% % % %         Y(loop,1:length(YinChunk))=YinChunk;
% % % %         freq=sampRate.*(1:L)/L; % frequency domain
% % % % 
% % % %         startIndex=startIndex+chunkSize-numOverlap;
% % % %         endIndex=endIndex+chunkSize-numOverlap;
% % % %         loop=loop+1;
% % % %     end
% % % %     endRampIndex=startIndex; %updates the endRampIndex to start where we ended the loop above
% % % %     meanY=mean(Y); %averages columns to find mean fft for the whole thing (averages the Ys of chunks)
% % % %     FzCorrFixed=ifft(meanY); %inverse fft to get back to amplitude of Fz
% % % %     plot(tDataCorr(startRampIndex:startRampIndex-1+length(FzCorrFixed),1,i),FzCorrFixed,'-k')
%     
%     startIndex=endRampIndex;
%     endIndex=startIndex+chunkSize;
%     loop=1;
%     Y=zeros(0);%ffts. each row is a chunk, each column a frequency.
%     while endIndex<(endRampIndex+40000) %loops through more data than will be graphed. dont need to do all of it
%         %to do: hann window
%         Fz=tDataCorr(startIndex:endIndex , 21, i); %the Fz in the chunk
%         L=length(Fz);%length of signal (in measurements)
%         YinChunk=fft(Fz,L); %fourier transform of the Fz
%         Y(loop,1:length(YinChunk))=YinChunk;
%         freq=sampRate.*(1:L)/L; % frequency domain
% 
%         startIndex=startIndex+chunkSize-numOverlap;
%         endIndex=endIndex+chunkSize-numOverlap;
%         loop=loop+1;
%     end
%     meanY=mean(Y); %averages columns to find mean fft for the whole thing (averages the Ys of chunks)
%     FzCorrFixed=ifft(meanY); %inverse fft to get back to amplitude of Fz
%     plot(tDataCorr(endRampIndex:endRampIndex-1+length(FzCorrFixed),1,i),FzCorrFixed,'-k')
    
    
    
    
    
%     %now trying to just do it manually... dont think I'll have better
%     %chances
%     totalLength=length(tDataCorr(:,21,i));
%     chunkSize=totalLength/1000;
%     loop=1;
%     startIndex=1;
%     endIndex=chunkSize;
%     PSD=zeros(0); %power spectral density. each row is one chunk, each column is a measurment
%     Y=zeros(0); %ffts. each row is a chunk, each column a freq
%     while endIndex<=totalLength %loops through all chunks
%         %maybe should apply a hann window?
%         Fz=tDataCorr(startIndex:endIndex , 21, i); %the Fz in the chunk
%         L=length(Fz);%length of signal (in measurements)
%         YinChunk=fft(Fz,L); %fourier transform of the Fz
%         Y(loop,1:length(YinChunk))=YinChunk;
%         PSDinChunk=YinChunk.*conj(YinChunk)/L;
%         PSD(loop,1:length(PSDinChunk))=PSDinChunk; %Power spectral density N^2/Hz (magnitude of the frequencies?)
%         freq=sampRate.*(1:L)/L; % frequency domain
%     
%         startIndex=startIndex+chunkSize;
%         endIndex=endIndex+chunkSize;
%         loop=loop+1;
%     end
%     meanY=mean(Y); %averages columns to find mean fft for the whole thing (averages the Ys of chunks)
%     FzCorrFixed=ifft(meanY); %inverse fft to get back to amplitude of Fz
%     plot(tDataCorr(1:chunkSize,1,i),FzCorrFixed,'-k')
%     %TAKEAWAY: this method (as I understand it, maybe I have the wrong idea entirely) does not work because the data does not fit
%     nicely into periodic functions. the data before ramp is very
%     different than the data during ramp is very different than data after
%     ramp. By breaking it into small chunks and averaging, I lose this.
%     SOLUTION: I guess I should chunk, do hanning window and fft on chunk,
%     make correction, inverse fft, then graph each chunk individually. 
%     POSSIBLY BETTER SOLUTION: Write code to determine time before ramp,
%     time during ramp, time after ramp. then can analyze each one of these reegions seperately.
%     by breaking into chunks, doing fft and averaging the ffts, then
%     inverse fft and plot.
    
    
    
    
% % % % % %     %trying to use pwelch method... I CANNOT FIGURE IT OJHT 
% % %     Fz=tDataCorr(:, 21, i);
% % %     wind=1000; %arbitrary-controlls how large they are.
% % %     %bandwidth=1; %??
% % %     %binwidth=2^nextpow2(sampRate/bandwidth);
% % %     %windowSize=hann(binwidth);
% % %     windowSize=length(Fz)/wind;
% % %     numOverlap=windowSize*0.5;
% % %     fre=1:sampRate; %IDK if this is correct. supposed to be every freq possible. aka all numbers of repetitions possible per unit time specified by sampRate
% % % % % %         %[Pxx, freqs]=pwelch(Fz, windowSize , numOverlap , fre , sampRate);
% % % %%%           [Pxx, freqs]=pwelch(Fz, windowSize , numOverlap , sampRate , 'onesided');
% % %     [Pxx,freqs]=pwelch(Fz,windowSize,numOverlap,[],sampRate);
% % %     fprintf('pwelched %f\n',i)
% % %     %amplitudeFFT contains real in column 1, imaginary in column 2
% % %     amplitudeFFT(1,:)=real(sqrt(Pxx*length(tDataCorr(length(tDataCorr(:,1,1)),1,i)))); %converts power spectral dnesity back to fft of amplitude (of Fz in this case) (by multiplying by T aka max time, sqrt rooting)
% % %     amplitudeFFT(2,:)=imag(sqrt(Pxx*length(tDataCorr(length(tDataCorr(:,1,1)),1,i))));
% % %     FzFixed=ifft(amplitudeFFT(1,:)+sqrt(amplitudeFFT(2,:).*conj(amplitudeFFT(2,:))))*sampRate; 
% % %     fprintf('%f\n%f\n',length(Pxx),length(amplitudeFFT))
% % %     %fprintf('%f\n%f',  , length(FzFixed))
% % %     %plot(tDataCorr(1+((chunk-1)*(totalLength)/1000) : chunk*totalLength/1000 ,1,i) , FzFixed,'-k')
% % %     plot(tDataCorr(1:length(FzFixed),1,i) , FzFixed,'-k') %just plotting at beginning arbitrarily. Later, will need to go back and plot it repeatedly from beginning to end.
% % %     
    
    
    
    
    %BELLOW: trying to make my own filter. kind of worked but not really.
%     L=length(tDataCorr(:,21,i));%length of signal (in measurements)
%     Y=fft(tDataCorr(:,21,i),L); %fourier transform of the Fz
%     PSD=Y.*conj(Y)/L; %Power spectral density N^2/Hz (magnitude of the frequencies?)
%     freq=sampRate.*(0:L)/L; % frequency domain
%     freq=freq(2:length(freq));
% %     plot(freq,PSD)
% %     title('PSD Stuff')
% %     xlabel('frequency (Hz)')
% %     ylabel('Power Spectral Density')
%     PSD999=prctile(PSD,99.9); %99.9th percentile of PSD
%     indices=PSD>PSD999; %array of 0s and 1s for if that index's PSD is in most important .01%. 1 is true, 0 false. MIGHT NEED TO ADD ANOTHER CONDITION/OTHER ARRAY FOR IF FREQ LARGE ENOUGH.
%     %PSDdenoised=PSD.*indices; %sets noise frequencies to have 0 power
%     Ydenoised=(Y.*indices); %sets noise frequencies to 0 in the fourier transform
%     FzCorrDenoised=ifft(Ydenoised); %inverse fft to get denoised signal
%     plot(tDataCorr(toDisplay,1,i),FzCorrDenoised(toDisplay),'-k')