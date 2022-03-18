%bestWindows is for the best windows, each row is a different test
%summedVals is an intermediary step. Each row is a different window, each
%column is a different test (after first 2 -starttime and endtime)
%windowsSorted is a sorted version of summedErrors. On top are the best
%windows (minimum sum). Each row is a different window. Columns are the 
%start/endtimes and the weighted sum. Each page is a different test.

function [bestWindows, summedVals, windowsSorted] = determineWindow(valueArray, weights)
    summedVals= [ valueArray(:,[1,2],1) , zeros(length(valueArray(:,1,1)),length(valueArray(1,1,:))) ]; %weighted sum of the error or std devs for the 7 metrics. each row is a different 10s window, each column is a different test (after first 2 rows, which are starttime and endtime)
    bestWindows=zeros(0); %columns are start and endtimes. Rows are different tests.
    for i=1:length(valueArray(1,1,:)) %loops through each test
        for k=1:length(valueArray(:,1,1)) %loops through each 10s window
            sum=0;
            for j=3:length(valueArray(1,:,1))%loops through each metric in 10s window
                sum=sum+(valueArray(k,j,i)*weights(j-2)); %sums up weighted error/st dev of different metrics
            end
            summedVals(k,i+2)=sum; %stores weighted sum of the error/st dev
        end
        [~,minIndex]=min(summedVals(:,i+2)); %finds index of lowest summed percent error
        bestWindows(i,[1,2])=valueArray(minIndex,[1,2],i); %uses this index to store the best 10s window
    end
    temporary=zeros(length(valueArray(:,1,1)),3,length(valueArray(1,1,:))); %layout is like windowsSorted, its just not sorted.
    for i=3:length(summedVals(1,:)) %loops through each test
        temporary(:,[1,2],i-2 )=summedVals(:,[1,2]); %adds starttimes and endtimes to temporary
        temporary(:,3,i-2)=summedVals(:,i); %adds everything else from the test
        windowsSorted(:,:,i-2)=sortrows(temporary(:,:,i-2),3); %sorts the test from smallest error on top, biggest on bottom
    end
    
    
end





