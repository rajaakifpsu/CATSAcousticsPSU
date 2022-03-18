clear
clc

%trial time 119 s
%0.00005 s between measurments

%Purpose: Analyze a single test from the research (tdhcs1) to see if I can
%come up with a method to extract frequency.
%I've saved the Full.txt to my computer, maybe later I can try to figure
%out how to import from onedrive so I dont need to download all of these...

M=readmatrix('tdhcs1.txt'); %importing data into matrix M

%METHOD 2 averages the time between when a blade is first picked up by the
%tacometer
soln=zeros(1,2);
for i=[4,5] %loops through 4th and 5th column of table, corresponding to motor1 and motor2 tacometer data
    measurmentsBetween=zeros(0);%array for number of measurements between starting and detecting first blade, between detections of new blades, and lastly between detection of last blade and end of measurments
    num=0; %number of measurments between new blades detected, recorded in mesurmentsBetween
    lastBlade=false; %boolean to describe if the tacometer sensed a blade in previous reading. (in the table, if previous reading was a 1)
    blade=false; %boolean to describe if the tacometer senses a blade in the current reading. 
    for j=1:height(M) %loops through every reading of the tacometer
        blade=M(j,i); %looks for if the blade is detected in this reading
        if (blade==true && lastBlade==false) %if it detects a new blade
            measurmentsBetween(length(measurmentsBetween)+1)=num; %records num
            num=0; %resets num to 0
        end
        lastBlade=blade;
        num=num+1;
    end
    avBtwBlade=mean(measurmentsBetween(2:(length(measurmentsBetween)))); %average number of measurments between detecting new blades (not including first measurment since that won't be consistent since it starts between blades, not on a blade)
    RPM=1/(avBtwBlade*0.00005/60); %converts to RPM, what we want
    %^Might need to add alittle if statement to prevent NaN results, make
    %them 0 instead
    soln(i-3)=RPM;
end

%results: it is 0 and 1894.7
%Should get: 0 and 1925. 30 off, but unsure if my code is wrong or if
%that's just the error involved...




% %METHOD 1 just counts the number of blades detected in the time the
% %tacometer is looking for them, then divides by time. I'm assuming the time
% %is 1s.
% fprintf('METHOD 1:\n')
% detected=0; %number of blades detected in the time the tacometer was looking. I'm assuming this to be 1 second.
% for i=[4,5] %loops through 4th and 5th columns of table, corresponding to motor1 and motor2 tacometer data
%     lastBlade=false; %boolean to describe if the tacometer sensed a blade in previous reading. (in the table, if previous reading was a 1)
%     blade=false; %boolean to describe if the tacometer senses a blade in the current reading. 
%     for j=1:height(M) %loops through every reading of the tacometer
%         blade=M(j,i); %looks for if the blade is detected in this reading
%         if (blade==true && lastBlade==false) %if it detects a new blade
%             detected=detected+1; %increment
%         end
%         lastBlade=blade;
%     end
%     detected=detected*60/119;
% end
