%% Bead plots final for paper

%load("'Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN13_16_STED12_BeadFitResults_ROI15_close8.mat")
fig17=figure(17);


plot(1:num2compare,amp_dymin16_mean/mean(amp_dymin16_mean(1:3)),'+g')
hold on
plot(1:num2compare,amp_dymin15_mean/mean(amp_dymin15_mean(1:3)),'ob')
plot(1:num2compare,amp_dymin14_mean/mean(amp_dymin14_mean(1:3)),'.c')
plot(1:num2compare,amp_dymin13_mean/mean(amp_dymin13_mean(1:3)),'xm')
plot(1:num2compare,amp_sted_mean/mean(amp_sted_mean(1:3)),'*r')


%Do poly fits to data for lines, order of 4 was what we agreed looked best
xx = [1:0.01:num2compare]';
p0 = polyfit(1:num2compare,amp_dymin16_mean/mean(amp_dymin16_mean(1:3)),4);
y0 = polyval(p0,xx);
plot(xx,y0,'-g')
p1 = polyfit(1:num2compare,amp_dymin15_mean/mean(amp_dymin15_mean(1:3)),4);
y1 = polyval(p1,xx);
plot(xx,y1,'-b')
p2 = polyfit(1:num2compare,amp_dymin14_mean/mean(amp_dymin14_mean(1:3)),4);
y2 = polyval(p2,xx);
plot(xx,y2,'-c')
p3 = polyfit(1:num2compare,amp_dymin13_mean/mean(amp_dymin13_mean(1:3)),4);
y3 = polyval(p3,xx);
plot(xx,y3,'-m')
p4 = polyfit(1:num2compare,amp_sted_mean/mean(amp_sted_mean(1:3)),4);
y4 = polyval(p4,xx);
plot(xx,y4,'-r')
xlabel("Image Number")
axis([1 num2compare 0 1.05])
ylabel("Mean Amplitude scaled at t=0")
legend( {'DyMIN1','DyMIN2','DyMIN3','DyMIN4','STED'}, 'Location',"best")
legend('boxoff')
fig17.Units               = 'centimeters';
fig17.Position(3)         = 8;
fig17.Position(4)         = 6;
set(fig17.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig17.PaperPositionMode   = 'auto';
print('Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\BeadMeanAmp_Poly', '-dpng', '-r600')
%% Plot num beads separately

numbeads_dymin16=numbeads_dymin16';
numbeads_dymin15=numbeads_dymin15';
numbeads_dymin14=numbeads_dymin14';
numbeads_dymin13=numbeads_dymin13';
numbeads_sted=numbeads_sted';
fig18=figure(18);
clf(fig18,'reset')
plot(1:num2compare,numbeads_dymin16/mean(numbeads_dymin16(1:3)),'+g')
hold on
plot(1:num2compare,numbeads_dymin15/mean(numbeads_dymin15(1:3)),'ob')
plot(1:num2compare,numbeads_dymin14/mean(numbeads_dymin14(1:3)),'.c')
plot(1:num2compare,numbeads_dymin13/mean(numbeads_dymin13(1:3)),'xm')
plot(1:num2compare,numbeads_sted/mean(numbeads_sted(1:3)),'*r')
% 
% 
% %Do poly fits to data for lines
xx = [1:0.01:num2compare]';
p0 = polyfit(1:num2compare,numbeads_dymin16/mean(numbeads_dymin16(1:3)),4);
y0 = polyval(p0,xx);
plot(xx,y0,'-g')
p1 = polyfit(1:num2compare,numbeads_dymin15/mean(numbeads_dymin15(1:3)),4);
y1 = polyval(p1,xx);
plot(xx,y1,'-b')
p2 = polyfit(1:num2compare,numbeads_dymin14/mean(numbeads_dymin14(1:3)),4);
y2 = polyval(p2,xx);
plot(xx,y2,'-c')
p3 = polyfit(1:num2compare,numbeads_dymin13/mean(numbeads_dymin13(1:3)),4);
y3 = polyval(p3,xx);
plot(xx,y3,'-m')
p4 = polyfit(1:num2compare,numbeads_sted/mean(numbeads_sted(1:3)),4);
y4 = polyval(p4,xx);
plot(xx,y4,'-r')
xlabel("Image Number")
ylabel("Number of beads scaled at t=0")
axis([1 num2compare 0 1.2])

fig18.Units               = 'centimeters';
fig18.Position(3)         = 8;
fig18.Position(4)         = 6;
set(fig18.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig18.PaperPositionMode   = 'auto';
print('Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\Numbeads_Poly', '-dpng', '-r600')

