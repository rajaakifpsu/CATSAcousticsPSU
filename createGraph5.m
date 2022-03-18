%in addition to graph 5 drawn in the file, I'll also graph the % Error by
%dividing error by the expected (the value in avFMRPMRep)

function createGraph5(errorsRep,FMRPMRep,fig)

% %graph as drawn in the file. - Don't care about. its wrong in file.
% figure(fig(1))
% hold on;
% grid on;
% 
% avErrors=mean(errorsRep,1); 
% maxErrors=max(errorsRep)-avErrors; 
% minErrors=avErrors-min(errorsRep); 
% 
% Xlabels = categorical({'Fx','Fy','Fz','Mx','My','Mz'});
% Xlabels = reordercats(Xlabels,{'Fx','Fy','Fz','Mx','My','Mz'});
% bar(Xlabels,avErrors(7:12),'r')
% er = errorbar(Xlabels,avErrors(7:12),minErrors(7:12),maxErrors(7:12)); %possibly incorrect?  
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% ylabel("Value Error (N or Nm)")
% title("Repeatability Error")
% hold off;



%new graph, percent error graph. for Fz Mz and RPM
figure(fig(2))
hold on;
grid on;
percErrors=zeros(size(errorsRep));
for i=1:length(percErrors(1,:))
    percErrors(:,i)=abs(errorsRep(:,i)./FMRPMRep(:,i))*100; %percent errors
end

avPercErrors=mean(percErrors,1); %average error for each metric (from the tests re-ran)
maxPercErrors=max(percErrors)-avPercErrors; %length of errorbar to max
minPercErrors=avPercErrors-min(percErrors); %length of errorbar to min

Xlabels = categorical({'Fz','Mz','RPM'});
Xlabels = reordercats(Xlabels,{'Fz','Mz','RPM'});
bar(Xlabels,avPercErrors([9,12,14]),'r')
er = errorbar(Xlabels,avPercErrors([9,12,14]),minPercErrors([9,12,14]),maxPercErrors([9,12,14]));    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ylabel("Percent Error")
ylim([0,20])
title("Repeatability Percent Error")
hold off;

end