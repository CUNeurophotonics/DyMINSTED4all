% First find beads (ROI) make sure that this is acceptable before moving on
% to the next step.
function [rois] = find_peaks(cdata, map, ROIspan, sigma, closefactor, edgebuf, warp_cut)

%% Find peaks and display their location
grayImage = cdata;
maxlim = max(max(cdata));
clims = [1 maxlim];%Specifies the limits for the min and max in the color plot of the image

p = FastPeakFind_bmh(grayImage, sigma); %edited to change sigma from main

num = size(p)/2.;
num = num(1);
%warp_cut fraction os frame size that should be cut off the left from warping
% warp_cut = .3;
k = 1;
for i=1:2:2*num
    valuesx(k)=p(i);
    valuesy(k)=p(i+1);
    k=k+1;
end

% figure(1);
% imagesc(cdata,clims)
% colormap(hot)
% hold on
% plot(valuesx,valuesy,'g+')
% set(h1, 'Position', [200 450 550 525])

%values_ are the original peaks
%val_ are the peaks that passed muster
valx = valuesx;
valy = valuesy;

%% Remove peaks that are too close together
%Selects those peaks that are too close to one another (under # pixels specified) and removes 

num = length(valuesx);
k=1;
test = 0;
for l=1:num
    for i=1:num
            if (l~=i) && (abs(sqrt((valuesy(l)-valuesy(i))^2+(valuesx(l)-valuesx(i))^2)) < closefactor)
                test(k)=i;
                test(k+1)=l;
                k=k+2;
            else
            end
    end
end
if test~=0
    test1 = unique(test);
    num = length(test1);
    k=0;
    for i = 1:num
        valx(test1(i)-k)=[];
        valy(test1(i)-k)=[];
        k=k+1;
    end
end

%% Gets rid of ROI's too near the edge of the frame
Imagesize = length(cdata);
k = 1;
for i=1:length(valx)
    x = valx(i);
    y = valy(i);
    if x<(edgebuf)
        badpointsx(k) = i;
        badpointsy(k) = i;
        k = k+1;
    elseif y<(edgebuf)
        badpointsx(k) = i;
        badpointsy(k) = i;
        k = k+1;
    elseif x >(Imagesize-(edgebuf))
        badpointsx(k) = i;
        badpointsy(k) = i;
        k = k+1;
    elseif y >(Imagesize-(edgebuf))
        badpointsx(k) = i;
        badpointsy(k) = i;
        k = k+1;
    elseif x < warp_cut*Imagesize
        badpointsx(k) = i;
        badpointsy(k) = i;
        k = k+1;
    elseif x > (1)*Imagesize  % was 1-warp_cut, but my images only have scanning weirdness on left hand side.
        badpointsx(k) = i;
        badpointsy(k) = i;
        k = k+1;
    end
end

k=0;
    if exist('badpointsx','var')  %ADDED 8/23/23 by Stephanie Pierce
        badpointsx
        for i = 1:length(badpointsx)
            valx(badpointsx(i)-k)=[];
            valy(badpointsy(i)-k)=[];
            k=k+1;
        end
    end
% Plots image and draws boxes after elimination
% h1 = figure(1)
% imagesc(cdata,clims)
% colormap(hot)
% hold on
% for i=1:length(valx)
%     rectangle('Position',[valx(i)-((ROIspan-1)/2+1),valy(i)-((ROIspan-1)/2+1),ROIspan,ROIspan],...
%         'LineWidth',2,'Edgecolor','g')
% end
% set(h1, 'Position', [800 450 350 350])
rois(1:length(valx), 1) = valx';
rois(1:length(valy), 2) = valy';
%keeps the points as ROI's if everything looks good
% button = questdlg('Keep these ROI`s?','ROI','Yes','Nope','Yes');
% switch button
%     case 'Yes'
%         rois(1:length(valx), 1) = valx';
%         rois(1:length(valy), 2) = valy';
%     case 'Nope'
%         disp('Try adjusting sigma.')
%         rois = 0;
end
