function createGraph10(errorsRep,FMRPMRep,errorsRepTotal,FMRPMRepTotal,fig)
loop=1;
for j=[1,4,7] %loops through the groups VDHCSAB BB CB
    percErrors=zeros(3,length(errorsRep(1,:))); %each row is a test, each column a metric
    for i=1:length(percErrors(1,:))
        percErrors(:,i)=abs(errorsRep(j:j+2,i)./FMRPMRep(j:j+2,i))*100; %percent errors
    end

    avPercErrors=mean(percErrors,1); %average error for each metric (from the tests re-ran)
    maxPercErrors=max(percErrors)-avPercErrors; %length of errorbar to max
    minPercErrors=avPercErrors-min(percErrors); %length of errorbar to min

    percErrorsTotal=zeros(3,length(errorsRepTotal(1,:)));
    for i=1:length(percErrorsTotal(1,:))
        percErrorsTotal(:,i)=abs(errorsRepTotal(j:j+2,i)./FMRPMRepTotal(j:j+2,i))*100;
    end

    avPercErrorsTot=mean(percErrorsTotal,1); %average error for each metric (from the tests re-ran)
    maxPercErrorsTot=max(percErrorsTotal)-avPercErrorsTot; %length of errorbar to max
    minPercErrorsTot=avPercErrorsTot-min(percErrorsTotal); %length of errorbar to min

    figure(fig(loop,1))%motor 1
    hold on;
    grid on;
    Xlabels = categorical({'Fz','Mz','RPM'});
    Xlabels = reordercats(Xlabels,{'Fz','Mz','RPM'});
    bar(Xlabels,avPercErrors([3,6,13]),'r')
    er = errorbar(Xlabels,avPercErrors([3,6,13]),minPercErrors([3,6,13]),maxPercErrors([3,6,13])); %bars for motor2 Fz Mz RPM 
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    ylabel("Percent Error")
    ylim([0,20])
    title("Repeatability Percent Error - Motor 1")
    hold off;

    figure(fig(loop, 2))%motor 2
    hold on;
    grid on;
    Xlabels = categorical({'Fz','Mz','RPM'});
    Xlabels = reordercats(Xlabels,{'Fz','Mz','RPM'});
    bar(Xlabels,avPercErrors([9,12,14]),'r')
    er = errorbar(Xlabels,avPercErrors([9,12,14]),minPercErrors([9,12,14]),maxPercErrors([9,12,14])); %bars for motor2 Fz Mz RPM 
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    ylabel("Percent Error")
    ylim([0,20])
    title("Repeatability Percent Error - Motor 2")
    hold off;

    figure(fig(loop, 3))%total
    hold on;
    grid on;
    Xlabels = categorical({'Fz','Mz','RPM'});
    Xlabels = reordercats(Xlabels,{'Fz','Mz','RPM'});
    bar(Xlabels,avPercErrorsTot([3,6,7]),'r')
    er = errorbar(Xlabels,avPercErrorsTot([3,6,7]),minPercErrorsTot([3,6,7]),maxPercErrorsTot([3,6,7])); %bars for motor2 Fz Mz RPM 
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    ylabel("Percent Error")
    ylim([0,20])
    title("Repeatability Percent Error - Total")
    hold off;
    
    loop=loop+1;
end
end