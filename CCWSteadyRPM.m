clear
clc

tests=["tdhcs1.txt","tdhcs2.txt","tdhcs3.txt","tdhcs4.txt","tdhcs5.txt","tdhcs6.txt","tdhcs7.txt","tdhcs15.txt","tdhcs16.txt"];
RPMs=zeros(length(tests),7);
%RPMs is matrix of solutions. first column is names - tdhcs___.
%2nd column is RPM of motor one.
%3rd column in RPM of motor two.
%Fourth column is desired RPM.
%5th column is error to desired
%6th column is measured RPM
%7th column is error to measured
RPMs(:,1)=[1:7,15,16]';
for i=1:length(tests)
   RPMs(i,[2,3])=getRPM(tests(i));
end
RPMs(:,4)=[1895,2493,2991,3490,3989,4487,4986,1496,1994]';
RPMs(:,5)=RPMs(:,4)-RPMs(:,3);
RPMs(:,6)=[1894,2484,3000,3477,4005,4463,4981,1510,2004]';
RPMs(:,7)=RPMs(:,6)-RPMs(:,3);

T=array2table(RPMs,'VariableNames',{'Name tdhcs__','m1 RPM', 'm2 RPM', 'desired RPM', 'to desired', 'measured RPM','to measured'})
writetable(T,'CCWSteadyRPMs.txt');

