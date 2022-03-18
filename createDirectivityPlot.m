%Purpose: to create a few directivity plot for oaspl comparing the unfiltered
%       data and the broadband data.

%inputs: (oaspl is vector representing overall sound pressure level for
%       unfiltered data for all mics
%       apl100 is that for sound pressure level in range [100, 20000]Hz 
%       spl500 is that for sound pressure level in range [500 20000]Hz
%       spl1k is the same but for range [1000 20000] Hz
%       spl10k is for range [10000 20000] Hz
%       )
%       ^For the unfiltered case, broadband (BB) and highpass only (HPO)
%       numofmics is the number of microphones
%       graphname is a string array of the name of the graphs to be made

%Outputs: a few graphs

%WARNING: NEED TO FIX GRAPHNAMES

function createDirectivityPlot(oaspl,oasplBB,oasplHPO,spl100,spl100BB,spl100HPO, spl500,spl500BB,spl500HPO,spl1k,spl1kBB,spl1kHPO,spl10k,spl10kBB,spl10kHPO, numofmics,graphname)

angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert into radians (only for plotting)

%Doing the directivity plots
% figure()
% polarplot(angle, oaspl) %the unfiltered plot
% thetalim([-40 , 30]);
% rlim([60 , 80]);

figure()%the unfiltered, broadband, and highpass only directivity plot for OASPL
polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'r',angle(1:(numofmics-1)),oasplBB(1:(numofmics-1)),'k',angle(1:(numofmics-1)),oasplHPO(1:(numofmics-1)),'--b'); %only first 11 since last mic is different size
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered','Broadband','Highpass Only')
title(graphname(1))

figure()%unfiltered, broadband, highpassonly directivity plot for spl in range [100, 20000]
polarplot(angle(1:(numofmics-1)),spl100(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl100BB(1:(numofmics-1)),'k',angle(1:(numofmics-1)),spl100HPO(1:(numofmics-1)),'--b'); %only first 11 since last mic is different size
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered','Broadband','Highpass Only')
title(graphname(2))

figure()%the unfiltered, broadband, and highpass only plot for spl in range [500, 20000]
polarplot(angle(1:(numofmics-1)),spl500(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl500BB(1:(numofmics-1)),'k',angle(1:(numofmics-1)),spl500HPO(1:(numofmics-1)),'--b')
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered','Broadband','Highpass Only')
title(graphname(3))

figure() %same for range [1000, 20000]
polarplot(angle(1:(numofmics-1)),spl1k(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl1kBB(1:(numofmics-1)),'k',angle(1:(numofmics-1)),spl1kHPO(1:(numofmics-1)),'--b')
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered','Broadband','Highpass Only')
title(graphname(4))

figure() %same for range [10000, 20000]
polarplot(angle(1:(numofmics-1)),spl10k(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl10kBB(1:(numofmics-1)),'k',angle(1:(numofmics-1)),spl10kHPO(1:(numofmics-1)),'--b')
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered','Broadband','Highpass Only')
title(graphname(5))

%M (plotting broadband in [100 20000] against unfiltered [0 20000 ]and highpassonly in [100 20000])
figure() 
polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl100HPO(1:(numofmics-1)),'--b',angle(1:(numofmics-1)),spl100BB(1:(numofmics-1)),'k')
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered OASPL','Highpass Only SPL [0.1, 20]kHz','Broadband SPL [0.1, 20]kHz')
title(graphname(6))

%N (BB [1000 20000], highpass [100 20000], unfiltered [0, 20000])
figure() %same but oaspl for unfiltered and highpassonly
polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl100HPO(1:(numofmics-1)),'--b',angle(1:(numofmics-1)),spl1kBB(1:(numofmics-1)),'k')
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered OASPL','Highpass Only SPL [0.1, 20]kHz','Broadband SPL [1, 20]kHz')
title(graphname(6))

figure() %O (BB [10000 20000], highpass [100 20000], unfiltered [0, 20000])
polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl100HPO(1:(numofmics-1)),'--b',angle(1:(numofmics-1)),spl10kBB(1:(numofmics-1)),'k')
thetalim([-40 , 30]);
rlim([60 , 80]);
legend('Unfiltered OASPL','Highpass Only SPL [0.1, 20]kHz','Broadband SPL [10, 20]kHz')
title(graphname(6))



%Bellow plot not useful
% figure()
% polarplot(angle(1:(numofmics-1)),oaspl(1:(numofmics-1)),'r',angle(1:(numofmics-1)),spl500BB(1:(numofmics-1)),'k')
% thetalim([-40 , 30]);
% rlim([60 , 80]);
% legend('Unfiltered OASPL','Broadband SPL range [0.5 20] kHz')
% title(graphname(3))





end