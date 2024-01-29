%% Plot FRC results  first load FRC data  "FRC Results All"
load('Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\FRC Results all.mat')
fig1 = figure(1);
plot(1:num2compare-1,1000*FRCResultsDyMIN16.FRCFixed17Calibratedmicron,'-+g')
hold on
plot(1:num2compare-1,1000*FRCResultsDyMIN15.FRCFixed17Calibratedmicron,'-ob')

plot(1:num2compare-1,1000*FRCResultsDyMIN14.FRCFixed17Calibratedmicron,'-.c')

plot(1:num2compare-1,1000*FRCResultsDyMIN13.FRCFixed17Calibratedmicron,'-xm')
plot(1:num2compare-1,1000*FRCResultsSTED12.FRCFixed17Calibratedmicron,'-*r')
d=2;

% Don't do polynomial fits, Emily thought it was better showing noise bc
% it's due to lower SNR
% xx = [1:0.01:num2compare-1]';
% x =[1:num2compare-1]';
% p0 = polyfit(x,1000*FRCResultsDyMIN16.FRCFixed17Calibratedmicron,d);
% y0 = polyval(p0,xx);
% plot(xx,y0,'-g')
% p1 = polyfit(x,1000*FRCResultsDyMIN15.FRCFixed17Calibratedmicron,d);
% y1 = polyval(p1,xx);
% plot(xx,y1,'-b')
% p2 = polyfit(x,1000*FRCResultsDyMIN14.FRCFixed17Calibratedmicron,d);
% y2 = polyval(p2,xx);
% plot(xx,y2,'-c')
% p3 = polyfit(x,1000*FRCResultsDyMIN13.FRCFixed17Calibratedmicron,d);
% y3 = polyval(p3,xx);
% plot(xx,y3,'-m')
% p4 = polyfit(x,1000*FRCResultsSTED12.FRCFixed17Calibratedmicron,d);
% y4 = polyval(p4,xx);
% plot(xx,y4,'-r')
xlabel("Image Number")
%axis([1 num2compare 0 1.05])
ylabel("Calibrated FRC 1/7")
% legend( {'DyMIN1','DyMIN2','DyMIN3','DyMIN4','STED'}, 'Location',"best")
% legend('boxoff')
axis([1 num2compare 50 120])

fig1.Units               = 'centimeters';
fig1.Position(3)         = 8;
fig1.Position(4)         = 6;
set(fig1.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig1.PaperPositionMode   = 'auto';
print('Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\FRCPlots', '-dpng', '-r600')
