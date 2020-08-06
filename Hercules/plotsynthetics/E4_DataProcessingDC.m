clear; close all; clc

load ./AuxFiles/EqkStnDB.mat

for nr = 1:size(Stn,2)
    load(['./FiguresSim2/Sim2DispCurves' num2str(nr) '.mat'])
    AllDCR(nr,:) = dC(1,:);
    AllDCT(nr,:) = dC(2,:);
    AllDCZ(nr,:) = dC(3,:);
end

AllDCR(abs(AllDCR) > 30) = 0;
AllDCT(abs(AllDCT) > 30) = 0;
AllDCZ(abs(AllDCZ) > 30) = 0;

subplot(1,3,1)
plot(Periods,mean(AllDCR),'.','linewidth',2), hold on
plot(Periods,mean(AllDCT),'.','linewidth',2), hold on
plot(Periods,mean(AllDCZ),'.','linewidth',2), hold on

xlim([4 52]), ylim([-3 3])
axis square; grid; xlabel('Period (s)'), ylabel('Mean \Deltac/c  (%)')
legend('R','T','Z')