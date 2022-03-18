%fig is 3xwhatever. each row is motor1,motor2,total. (whatever=number of
%tests)
%groupLen is number of tests in a group of tests.
%xlimit. each row is motor1 for tests then motor2 for tests, then total for
%tests (all tests within 1 group)
%, each column a limit min then max
%, each page is a different group
%chunkPortion, overlap used in createGraph4
%titles used for titles duh. each row is a different test, columns are for
%motor 1 motor2 and total

function createGraph9(RPMs1,RPMs2,tData,tDataCorr, xlimit, sampRate, groupLen, chunkPortion, overlap, titles, fig)
loop=0;
for i=1:groupLen:length(tDataCorr(1,1,:)) %loops through all the groups

    %NOTE: when passing parameters to createGraph4, will need to be careful
    %with fig. Because I coded createGraph4 stupidly, it uses the first numbers
    %specified by fig to plot motor2 stuff first, then itll go to motor1. To
    %counteract how this is reversed, I'll pass the fig parameter in reversed
    %order, so hopefully itll come out right.

    %call createGraph4 on the full data, do trick with fig
    passFig=fig([2,1],i:i+groupLen-1);
    passXLimit=xlimit(1:(2*groupLen) , :,loop+1); %x limits for motor1 and 2, in the current group
    passTitles=titles(i:i+groupLen-1,[1,2]);
    createGraph4(RPMs1(:,:,i:i+groupLen-1),RPMs2(:,:,i:i+groupLen-1),tData(:,:,i:i+groupLen-1),tDataCorr(:,:,i:i+groupLen-1), passXLimit, sampRate,groupLen,chunkPortion, overlap, passTitles, passFig);

    %create new data for combined case. Add the forces and moments. %make RPMs1
    %empty so createGraph4 wont create unneeded plots. Simply use RPMs2 there.
    %Then, in this funciton, graph RPMs1 on top of it. 

    totalTDataCorr=zeros(size(tDataCorr));
    totalTDataCorr(:,1,:)=tDataCorr(:,1,:); %the times is the same.
    totalTDataCorr(:,[21,24],:)=tDataCorr(:,[21,24],:)+tDataCorr(:,[15,18],:); %gets the total Fz and Mz columns
    RPMs1tot=[];
    RPMs2tot=zeros(2,2,length(RPMs2(1,1,:)));%RPMs2;
    passFig=fig(3,i:i+groupLen-1);
    L=length(xlimit(:,1));
    passXLimit=xlimit(L-groupLen+1:L,:,loop+1); %limits for total case in the current group
    passTitles=[repmat("",groupLen,1) , titles(i:i+groupLen-1,3)]; %needs to be in 2nd col since the createGraph4 function will interpret it if it were motor2
    createGraph4(RPMs1tot,RPMs2tot,totalTDataCorr(:,:,i:i+groupLen-1),totalTDataCorr(:,:,i:i+groupLen-1),passXLimit,sampRate,groupLen,chunkPortion, overlap, passTitles, passFig);
    
    for j=1:groupLen %loops through all tests in the group
        figure (fig(3,j+(loop*groupLen))) %selects the figure it for that test, for the total graph
        subplot(3,1,3) %the rpm plot
        hold on;
        stop=find( RPMs1(1,:,(j+loop*groupLen))==0 , 1,'first')-1; %looks for when RPM2 goes to zero in the matrix and stops graphing after that point (since thats when the data stopped being collected. after that is just a bunch of 0s that need to be there since each test will have a different number of RPMs recorded)
        if isempty(stop)
            stop=length(RPMs1(1,:,j+(loop*groupLen))); %if the RPM never goes to zero (has recordings for all positions in the matrix), then stop at the end
        end
        plot(RPMs1(2,1:stop,j+(loop*groupLen)) , RPMs1(1,1:stop,j+(loop*groupLen)),'-k')
        
        stop=find( RPMs2(1,:,j+(loop*groupLen))==0 , 1,'first')-1; %looks for when RPM2 goes to zero in the matrix and stops graphing after that point (since thats when the data stopped being collected. after that is just a bunch of 0s that need to be there since each test will have a different number of RPMs recorded)
        if isempty(stop)
            stop=length(RPMs2(1,:,j+(loop*groupLen))); %if the RPM never goes to zero (has recordings for all positions in the matrix), then stop at the end
        end
        plot(RPMs2(2,1:stop,j+(loop*groupLen)) , RPMs2(1,1:stop,j+(loop*groupLen)),'-r')
        legend('Motor 1','Motor 2','Location','northeast')
        hold off;
        
    end
    
    loop=loop+1;
end

end