%COMMENT TEST

function createGraphs7n8(FMRPM,fig)

S=["No difference","Low difference","High difference"];
for i=1:3 %getting strings for labels
    S(i)=strcat(S(i)," (",string(round(FMRPM(i,13)))," & ",string(round(FMRPM(i,14)))," RPM)" );
end    

X=categorical(S);

grid on;
figure(fig(1))
bar(X,FMRPM(:,[3,9]))
ylabel('Thrust (N)')
title("Thrust")
legend("Motor1","Motor2","Location","northeastoutside")

figure(fig(2))
grid on;
bar(X,FMRPM(:,[6,12]))
ylabel('Torque (Nm)')
title("Torque")
legend("Motor1","Motor2","Location","northeastoutside")

end