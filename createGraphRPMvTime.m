%singleTestData is the test data from a single test. ie one page, not the
%whole 3d thing

function createGraphRPMvTime(singleTestData, sampRate,fig)

figure(fig)
[~, RPMmot2]=getRPMRamp(singleTestData, sampRate); %not ramp, but this function is useful regardless
scatter(RPMmot2(2,:),RPMmot2(1,:),'kx')
ylabel("RPM")
xlabel("Time (s)")
title("RPM v Time (T-Motor ~5000 RPM)")

end