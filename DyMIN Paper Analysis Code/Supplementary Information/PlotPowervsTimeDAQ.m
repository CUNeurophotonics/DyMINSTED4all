%%%%  Find trigger threshold, start t=0 there, create table, convert
%%%%  duration to seconds
load('DAQ_10mW_Pvst.mat')
calibrate = 59.2613;
temp = find(DAQ_10mW.ai3>1.5 );
temp = temp(1,1)-3;

DAQ_10mW_t = timetable2table(DAQ_10mW,'ConvertRowTimes',false);
DAQ_10mW_t = removevars(DAQ_10mW_t,'ai3');
DAQ_10mW_t.Time = seconds(DAQ_10mW.Time);
%Tnew([18,20,21],:) = [];
DAQ_10mW_t(1:temp,:)=[];
DAQ_10mW_t.Time =DAQ_10mW_t.Time-DAQ_10mW_t.Time(1);


%% The data is way too noisy, average every 10 data points from 10^-4 s and up
temp0 = find(DAQ_10mW_t.Time>10^-4);
temp0 = temp0(1,1);
temp1 = find(DAQ_10mW_t.Time>10^-2);
temp1 = temp1(1,1);
temp2 = length(DAQ_10mW_t.ai2);
data0 = DAQ_10mW_t.ai2(1:temp0);
data1 = medfilt1(DAQ_10mW_t.ai2(temp0:temp1+1),10);
%y = medfilt1(x,n)
data2 = medfilt1(DAQ_10mW_t.ai2(temp1:temp2),100);
DAQ_10mW_t.medfilt = cat(1,data0,data1(2:length(data1)-1),data2(2:length(data2))).*calibrate;
DAQ_10mW_f = removevars(DAQ_10mW_t,'ai2');
save("DAQ_10mW_Pvst_filt_cal.mat", "DAQ_10mW_f")
%% repeat for 20mW data
load('DAQ_20mW_Pvst.mat')
temp = find(DAQ_20mW.ai3>1.5);
temp = temp(1,1)-3;

DAQ_20mW_t = timetable2table(DAQ_20mW,'ConvertRowTimes',false);
DAQ_20mW_t = removevars(DAQ_20mW_t,'ai3');
DAQ_20mW_t.Time = seconds(DAQ_20mW.Time);
%Tnew([18,20,21],:) = [];
DAQ_20mW_t(1:temp,:)=[];
DAQ_20mW_t.Time =DAQ_20mW_t.Time-DAQ_20mW_t.Time(1);

temp0 = find(DAQ_20mW_t.Time>10^-4);
temp0 = temp0(1,1);
temp1 = find(DAQ_20mW_t.Time>10^-2);
temp1 = temp1(1,1);
temp2 = length(DAQ_20mW_t.ai2);
data0 = DAQ_20mW_t.ai2(1:temp0);
data1 = medfilt1(DAQ_20mW_t.ai2(temp0:temp1+1),10);
%y = medfilt1(x,n)
data2 = medfilt1(DAQ_20mW_t.ai2(temp1:temp2),100);
DAQ_20mW_t.medfilt = cat(1,data0,data1(2:length(data1)-1),data2(2:length(data2))).*calibrate;
DAQ_20mW_f = removevars(DAQ_20mW_t,'ai2');
save("DAQ_20mW_Pvst_filt_cal.mat", "DAQ_20mW_f")
%% repeat for 30mW data
load('DAQ_30mW_Pvst.mat')
temp = find(DAQ_30mW.ai3>1.5);
temp = temp(1,1)-3;

DAQ_30mW_t = timetable2table(DAQ_30mW,'ConvertRowTimes',false);
DAQ_30mW_t = removevars(DAQ_30mW_t,'ai3');
DAQ_30mW_t.Time = seconds(DAQ_30mW.Time);
%Tnew([18,20,21],:) = [];
DAQ_30mW_t(1:temp,:)=[];
DAQ_30mW_t.Time =DAQ_30mW_t.Time-DAQ_30mW_t.Time(1);

temp0 = find(DAQ_30mW_t.Time>10^-4);
temp0 = temp0(1,1);
temp1 = find(DAQ_30mW_t.Time>10^-2);
temp1 = temp1(1,1);
temp2 = length(DAQ_30mW_t.ai2);
data0 = DAQ_30mW_t.ai2(1:temp0);
data1 = medfilt1(DAQ_30mW_t.ai2(temp0:temp1+1),10);
%y = medfilt1(x,n)
data2 = medfilt1(DAQ_30mW_t.ai2(temp1:temp2),100);
DAQ_30mW_t.medfilt = cat(1,data0,data1(2:length(data1)-1),data2(2:length(data2))).*calibrate;
DAQ_30mW_f = removevars(DAQ_30mW_t,'ai2');
save("DAQ_30mW_Pvst_filt_cal.mat", "DAQ_30mW_f")

%% repeat for 40mW data
load('DAQ_40mW_Pvst.mat')
temp = find(DAQ_40mW.ai3>1.5);
temp = temp(1,1)-3;

DAQ_40mW_t = timetable2table(DAQ_40mW,'ConvertRowTimes',false);
DAQ_40mW_t = removevars(DAQ_40mW_t,'ai3');
DAQ_40mW_t.Time = seconds(DAQ_40mW.Time);
%Tnew([18,20,21],:) = [];
DAQ_40mW_t(1:temp,:)=[];
DAQ_40mW_t.Time =DAQ_40mW_t.Time-DAQ_40mW_t.Time(1);

temp0 = find(DAQ_40mW_t.Time>10^-4);
temp0 = temp0(1,1);
temp1 = find(DAQ_40mW_t.Time>10^-2);
temp1 = temp1(1,1);
temp2 = length(DAQ_40mW_t.ai2);
data0 = DAQ_40mW_t.ai2(1:temp0);
data1 = medfilt1(DAQ_40mW_t.ai2(temp0:temp1+1),10);
%y = medfilt1(x,n)
data2 = medfilt1(DAQ_40mW_t.ai2(temp1:temp2),100);
DAQ_40mW_t.medfilt = cat(1,data0,data1(2:length(data1)-1),data2(2:length(data2))).*calibrate;
DAQ_40mW_f = removevars(DAQ_40mW_t,'ai2');
save("DAQ_40mW_Pvst_filt_cal.mat", "DAQ_40mW_f")
%% repeat for 50mW data
load('DAQ_50mW_Pvst.mat')
temp = find(DAQ_50mW.ai3>1.5);
temp = temp(1,1)-3;

DAQ_50mW_t = timetable2table(DAQ_50mW,'ConvertRowTimes',false);
DAQ_50mW_t = removevars(DAQ_50mW_t,'ai3');
DAQ_50mW_t.Time = seconds(DAQ_50mW.Time);
%Tnew([18,20,21],:) = [];
DAQ_50mW_t(1:temp,:)=[];
DAQ_50mW_t.Time =DAQ_50mW_t.Time-DAQ_50mW_t.Time(1);

temp0 = find(DAQ_50mW_t.Time>10^-4);
temp0 = temp0(1,1);
temp1 = find(DAQ_50mW_t.Time>10^-2);
temp1 = temp1(1,1);
temp2 = length(DAQ_50mW_t.ai2);
data0 = DAQ_50mW_t.ai2(1:temp0);
data1 = medfilt1(DAQ_50mW_t.ai2(temp0:temp1+1),10);
%y = medfilt1(x,n)
data2 = medfilt1(DAQ_50mW_t.ai2(temp1:temp2),100);
DAQ_50mW_t.medfilt = cat(1,data0,data1(2:length(data1)-1),data2(2:length(data2))).*calibrate;
DAQ_50mW_f = removevars(DAQ_50mW_t,'ai2');
save("DAQ_50mW_Pvst_filt_cal.mat", "DAQ_50mW_f")
%% repeat for 60mW data
load('DAQ_60mW_Pvst.mat')
temp = find(DAQ_60mW.ai3>1.5);
temp = temp(1,1)-3;

DAQ_60mW_t = timetable2table(DAQ_60mW,'ConvertRowTimes',false);
DAQ_60mW_t = removevars(DAQ_60mW_t,'ai3');
DAQ_60mW_t.Time = seconds(DAQ_60mW.Time);
%Tnew([18,20,21],:) = [];
DAQ_60mW_t(1:temp,:)=[];
DAQ_60mW_t.Time =DAQ_60mW_t.Time-DAQ_60mW_t.Time(1);

temp0 = find(DAQ_60mW_t.Time>10^-4);
temp0 = temp0(1,1);
temp1 = find(DAQ_60mW_t.Time>10^-2);
temp1 = temp1(1,1);
temp2 = length(DAQ_60mW_t.ai2);
data0 = DAQ_60mW_t.ai2(1:temp0);
data1 = medfilt1(DAQ_60mW_t.ai2(temp0:temp1+1),10);
%y = medfilt1(x,n)
data2 = medfilt1(DAQ_60mW_t.ai2(temp1:temp2),100);
DAQ_60mW_t.medfilt = cat(1,data0,data1(2:length(data1)-1),data2(2:length(data2))).*calibrate;
DAQ_60mW_f = removevars(DAQ_60mW_t,'ai2');
save("DAQ_60mW_Pvst_filt_cal.mat", "DAQ_60mW_f")
%% Clear timetables
clear DAQ_10mW DAQ_20mW DAQ_30mW DAQ_40mW DAQ_50mW DAQ_60mW

%% 
figure(1)
subplot(1,2,1)
semilogx(DAQ_10mW_f.Time,DAQ_10mW_f.medfilt)
hold on
semilogx(DAQ_20mW_f.Time,DAQ_20mW_f.medfilt)
semilogx(DAQ_30mW_f.Time,DAQ_30mW_f.medfilt)
semilogx(DAQ_40mW_f.Time,DAQ_40mW_f.medfilt)
semilogx(DAQ_50mW_f.Time,DAQ_50mW_f.medfilt)
semilogx(DAQ_60mW_f.Time,DAQ_60mW_f.medfilt)
hold off
xlabel('Time (s)')
ylabel('Power (mW)')
axis([10^-6 10^2 0 70])
%% Now calculate the final power for the 60mW trace and find the corresponding PSTED effective.  
% Then I can use thee effective STED power for the comparison conventional
% STED images, but set the Pfinal for the STED power during DyMIN.  Do the
% same for whatever power I decide on for P2 that corresponds to that
% STED value.
P3effective = zeros(6,1);
P3final = zeros(6,1);
index1 = find(DAQ_60mW_f.Time>=10^-5,1);
index2 = find(DAQ_60mW_f.Time>=10^-4,1);
P3effective(6,1) = mean(DAQ_60mW_f.medfilt(index1:index2));
index3 = find(DAQ_60mW_f.Time>=60,1);
index4 = size(DAQ_60mW_f.Time,1);
P3final(6,1) = mean(DAQ_60mW_f.medfilt(index3:index4));
% Do the same for every power value and see if I can get a linear fit
% I need to calculate Pset for step 2 to get 4mW effective power
index1 = find(DAQ_50mW_f.Time>=10^-5,1);
index2 = find(DAQ_50mW_f.Time>=10^-4,1);
P3effective(5,1) = mean(DAQ_50mW_f.medfilt(index1:index2));
index3 = find(DAQ_50mW_f.Time>=60,1);
index4 = size(DAQ_50mW_f.Time,1);
P3final(5,1) = mean(DAQ_50mW_f.medfilt(index3:index4));

index1 = find(DAQ_40mW_f.Time>=10^-5,1);
index2 = find(DAQ_40mW_f.Time>=10^-4,1);
P3effective(4,1) = mean(DAQ_40mW_f.medfilt(index1:index2));
index3 = find(DAQ_40mW_f.Time>=60,1);
index4 = size(DAQ_40mW_f.Time,1);
P3final(4,1) = mean(DAQ_40mW_f.medfilt(index3:index4));

index1 = find(DAQ_30mW_f.Time>=10^-5,1);
index2 = find(DAQ_30mW_f.Time>=10^-4,1);
P3effective(3,1) = mean(DAQ_30mW_f.medfilt(index1:index2));
index3 = find(DAQ_30mW_f.Time>=60,1);
index4 = size(DAQ_30mW_f.Time,1);
P3final(3,1) = mean(DAQ_30mW_f.medfilt(index3:index4));

index1 = find(DAQ_20mW_f.Time>=10^-5,1);
index2 = find(DAQ_20mW_f.Time>=10^-4,1);
P3effective(2,1) = mean(DAQ_20mW_f.medfilt(index1:index2));
index3 = find(DAQ_20mW_f.Time>=60,1);
index4 = size(DAQ_20mW_f.Time,1);
P3final(2,1) = mean(DAQ_20mW_f.medfilt(index3:index4));

index1 = find(DAQ_10mW_f.Time>=10^-5,1);
index2 = find(DAQ_10mW_f.Time>=10^-4,1);
P3effective(1,1) = mean(DAQ_10mW_f.medfilt(index1:index2));
index3 = find(DAQ_10mW_f.Time>=60,1);
index4 = size(DAQ_10mW_f.Time,1);
P3final(1,1) = mean(DAQ_10mW_f.medfilt(index3:index4));
%% Plotting and fitting

subplot(1,2,2)

f=fit(P3effective,P3final,'poly2')
plot(f,P3effective, P3final)
xlabel('Initial Power (mW)')
ylabel('Final Power (mW)')

