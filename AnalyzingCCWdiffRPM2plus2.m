clearvars -except testData
clc
close all

tests=["vdhcsab7","vdhclab7","vdhchab7","vdhcsbb7","vdhclbb7","vdhchbb7","vdhcscb7","vdhclcb7","vdhchcb7"]; %in groups of three by seperation distance. steady, low difference, high difference (in RPM)
testnums=7:.1:7.8;
sampleRate=20000; %20000 measurements/s

%imports all the data from all the txt files into a 3-dimensional array.
%columns and rows correspond to a single txt file, different pages (3rd index)are
%different txt files.
if ~exist('testData')
    testData=zeros(580000,27,length(tests));
    for i=1:length(tests)
        testData(:,:,i)=readmatrix(tests(i));
        fprintf('read %f\n',i)
    end
    save('testDataVdhcDiffRPM','testData','-v7.3');
else
    fprintf('File data already stored\n')
end



