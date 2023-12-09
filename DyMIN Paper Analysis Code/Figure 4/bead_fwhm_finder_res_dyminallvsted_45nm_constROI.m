% Cut together all the image analysis stuff Steph gave me into one thing
% First find beads (ROI) make sure that this is acceptable before moving on
% to the next step. Then fit all the beads to gaussians and make sure each
% of these is acceptable. Loop through this routine for multiple images


% if noise is being read as a peak, adjust the threshold in
% 'FastPeakFind_bmh'
clearvars
close all
%% Define path, Add folders with images and generate list
%add path to code and path to data folders
addpath(genpath('Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED_Matlab\'))
addpath("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN13\");
addpath("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN14\");
addpath("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN15\");
addpath("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN16\");
addpath("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\STED12\");
%naming information on image to analyze
STEDfilelist = dir("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\STED12");
DyMIN13filelist = dir("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN13");
DyMIN14filelist = dir("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN14");
DyMIN15filelist = dir("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN15");
DyMIN16filelist = dir("Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\DyMIN16");
%remove first two lines that are just invisible files no data
STEDfilelist(1:3) = [];
DyMIN13filelist(1:2) = [];
DyMIN14filelist(1:2) = [];
DyMIN15filelist(1:2) = [];
DyMIN16filelist(1:2) = [];
STEDfilelist = {STEDfilelist.name}.';
DyMIN13filelist = {DyMIN13filelist.name}.';
DyMIN14filelist = {DyMIN14filelist.name}.';
DyMIN15filelist = {DyMIN15filelist.name}.';
DyMIN16filelist = {DyMIN16filelist.name}.';
[num2compare, blah] = size(STEDfilelist);  %num2compare is the number of timelapse images, or number of power data points for res vs P measurements

%% More parameters

%These are cell arrays, index them using STEDfilelist{1,1} through
%STEDfilelist{20,1}

%parameter of image
frame_length = 10; %length of one side of image in microns
Npixels = 512; % number of pixel for one side of image
maxnumbeads = 215;  %In order to generate arrays to save, assume I have <=this many beads in image
nroi2compare=1;     %Added by SP 8/24/23, so I can easily pass only 1 roi and use this same code for comparing dymin vs sted
%% parmeters for peak finding
sigma = 1.8; %initial guess at FWHM, in pixels

%num of pix on each side of ROI center pixel, 
%total length is of ROI is 2*pixel+1
ROIspanmax = 15; %needs to be odd, max # of pixels for ROI span
closefactor = 8; % pixel spacing from one ROI to another
edgebuf = (ROIspanmax-1)/2+1; % needs to be at least 1 pixel
rois_sted = zeros(maxnumbeads, 2, num2compare);
rois_dymin13 = zeros(maxnumbeads, 2, num2compare);
rois_dymin14 = zeros(maxnumbeads, 2, num2compare);
rois_dymin15 = zeros(maxnumbeads, 2, num2compare);
rois_dymin16 = zeros(maxnumbeads, 2, num2compare);
pixelsize = frame_length/Npixels;  %in microns
mingoodness = .6;
warp_cut = 0.1; % don't use the left or right warp_cut% of the image due to scan issues

% these parameters can be used to get rid of fits smaller or larger than is
% reasonable (i.e. can't resolve something smaller than pixel and don't
FWHMmax = .4; %throw out spots larger than about 400 nm
FWHMmin = .02; %throw out fits smaller than 20 nm

%% Initialize the arrays to hold the fit parameters and results in
choice = 1; %this is a dialog choice of whether to accept the ROIs
fwhm_x_median_sted = zeros(1, num2compare);
fwhm_y_median_sted = zeros(1, num2compare);
numbeads_sted = zeros(1, num2compare);
fwhm_x_all_sted = zeros(maxnumbeads,num2compare);
fwhm_y_all_sted = zeros(maxnumbeads,num2compare);
fwhm_xyav_all_sted = zeros(maxnumbeads,num2compare);
fwhm_unfilt_sted = zeros(maxnumbeads,2,num2compare);
ROIgoodness_all_sted = zeros(2,nroi2compare,num2compare);
fwhm_unfilt_sted_all = zeros(maxnumbeads,2,20);
median_sted_all = zeros(num2compare,3);
mean_sted_all = zeros(num2compare,3);
std_sted_all = zeros(num2compare,3);
amp_sted_all = zeros(maxnumbeads,num2compare);
bg_sted_all = zeros(maxnumbeads,num2compare);

%% Initialize variables for each DyMIN trace I am going to compare

fwhm_x_median_dymin13 = zeros(1, num2compare);
fwhm_y_median_dymin13 = zeros(1, num2compare);
numbeads_dymin13 = zeros(1, num2compare);
fwhm_x_all_dymin13 = zeros(maxnumbeads,num2compare);
fwhm_y_all_dymin13 = zeros(maxnumbeads,num2compare);
fwhm_xyav_all_dymin13 = zeros(maxnumbeads,num2compare);
fwhm_unfilt_dymin13 = zeros(maxnumbeads,2,num2compare);
ROIgoodness_all_dymin13 = zeros(2,nroi2compare,num2compare);
fwhm_unfilt_dymin13_all = zeros(maxnumbeads,2,20);
median_dymin13_all = zeros(num2compare,3);
mean_dymin13_all = zeros(num2compare,3);
std_dymin13_all = zeros(num2compare,3);
amp_dymin13_all = zeros(maxnumbeads,num2compare);
bg_dymin13_all = zeros(maxnumbeads,num2compare);
STEDhalfROIused = zeros(1,num2compare);
DyMIN13halfROIused =zeros(1,num2compare);

fwhm_x_median_dymin14 = zeros(1, num2compare);
fwhm_y_median_dymin14 = zeros(1, num2compare);
numbeads_dymin14 = zeros(1, num2compare);
fwhm_x_all_dymin14 = zeros(maxnumbeads,num2compare);
fwhm_y_all_dymin14 = zeros(maxnumbeads,num2compare);
fwhm_xyav_all_dymin14 = zeros(maxnumbeads,num2compare);
fwhm_unfilt_dymin14 = zeros(maxnumbeads,2,num2compare);
ROIgoodness_all_dymin14 = zeros(2,nroi2compare,num2compare);
fwhm_unfilt_dymin14_all = zeros(maxnumbeads,2,20);
median_dymin14_all = zeros(num2compare,3);
mean_dymin14_all = zeros(num2compare,3);
std_dymin14_all = zeros(num2compare,3);
amp_dymin14_all = zeros(maxnumbeads,num2compare);
bg_dymin14_all = zeros(maxnumbeads,num2compare);
DyMIN14halfROIused =zeros(1,num2compare);
fwhm_x_median_dymin14 = zeros(1, num2compare);
fwhm_y_median_dymin14 = zeros(1, num2compare);
numbeads_dymin14 = zeros(1, num2compare);
fwhm_x_all_dymin14 = zeros(maxnumbeads,num2compare);
fwhm_y_all_dymin14 = zeros(maxnumbeads,num2compare);
fwhm_xyav_all_dymin14 = zeros(maxnumbeads,num2compare);
fwhm_unfilt_dymin14 = zeros(maxnumbeads,2,num2compare);
ROIgoodness_all_dymin14 = zeros(2,nroi2compare,num2compare);
fwhm_unfilt_dymin14_all = zeros(maxnumbeads,2,20);
median_dymin14_all = zeros(num2compare,3);
mean_dymin14_all = zeros(num2compare,3);
std_dymin14_all = zeros(num2compare,3);
amp_dymin14_all = zeros(maxnumbeads,num2compare);
bg_dymin14_all = zeros(maxnumbeads,num2compare);
DyMIN14halfROIused =zeros(1,num2compare);

fwhm_x_median_dymin15 = zeros(1, num2compare);
fwhm_y_median_dymin15 = zeros(1, num2compare);
numbeads_dymin15 = zeros(1, num2compare);
fwhm_x_all_dymin15 = zeros(maxnumbeads,num2compare);
fwhm_y_all_dymin15 = zeros(maxnumbeads,num2compare);
fwhm_xyav_all_dymin15 = zeros(maxnumbeads,num2compare);
fwhm_unfilt_dymin15 = zeros(maxnumbeads,2,num2compare);
ROIgoodness_all_dymin15 = zeros(2,nroi2compare,num2compare);
fwhm_unfilt_dymin15_all = zeros(maxnumbeads,2,20);
median_dymin15_all = zeros(num2compare,3);
mean_dymin15_all = zeros(num2compare,3);
std_dymin15_all = zeros(num2compare,3);
amp_dymin15_all = zeros(maxnumbeads,num2compare);
bg_dymin15_all = zeros(maxnumbeads,num2compare);
DyMIN15halfROIused =zeros(1,num2compare);

fwhm_x_median_dymin16 = zeros(1, num2compare);
fwhm_y_median_dymin16 = zeros(1, num2compare);
numbeads_dymin16 = zeros(1, num2compare);
fwhm_x_all_dymin16 = zeros(maxnumbeads,num2compare);
fwhm_y_all_dymin16 = zeros(maxnumbeads,num2compare);
fwhm_xyav_all_dymin16 = zeros(maxnumbeads,num2compare);
fwhm_unfilt_dymin16 = zeros(maxnumbeads,2,num2compare);
ROIgoodness_all_dymin16 = zeros(2,nroi2compare,num2compare);
fwhm_unfilt_dymin16_all = zeros(maxnumbeads,2,20);
median_dymin16_all = zeros(num2compare,3);
mean_dymin16_all = zeros(num2compare,3);
std_dymin16_all = zeros(num2compare,3);
amp_dymin16_all = zeros(maxnumbeads,num2compare);
bg_dymin16_all = zeros(maxnumbeads,num2compare);
DyMIN16halfROIused =zeros(1,num2compare);

tag = 'sted';

%% loop through, load images, find peaks, fit to Lorentzian
for ii=1:num2compare
    [cdata_sted,map_sted] = imread(STEDfilelist{ii,1});
    [cdata_dymin13,map_dymin13] = imread(DyMIN13filelist{ii,1});
    [cdata_dymin14,map_dymin14] = imread(DyMIN14filelist{ii,1});
    [cdata_dymin15,map_dymin15] = imread(DyMIN15filelist{ii,1});
    [cdata_dymin16,map_dymin16] = imread(DyMIN16filelist{ii,1});
    %find peaks of image (regions of interest) in conventional sted image
    [rois_sted] = find_peaks(cdata_sted, map_sted, ROIspanmax, sigma, closefactor, edgebuf, warp_cut);
    [rois_dymin13] = find_peaks(cdata_dymin13, map_dymin13, ROIspanmax, sigma, closefactor, edgebuf, warp_cut);
    [rois_dymin14] = find_peaks(cdata_dymin14, map_dymin14, ROIspanmax, sigma, closefactor, edgebuf, warp_cut);
    [rois_dymin15] = find_peaks(cdata_dymin15, map_dymin15, ROIspanmax, sigma, closefactor, edgebuf, warp_cut);
    [rois_dymin16] = find_peaks(cdata_dymin16, map_dymin16, ROIspanmax, sigma, closefactor, edgebuf, warp_cut);
    if rois_sted ~= 0
        %[fwhm_x, fwhm_y,fwhm_unfilt, amp, bg, ROIgoodnessinfo] =fwhm_fit(cdata, ...
        % ~, ROIspan, pixelsize, rois, mingoodness, FWHMmax, FWHMmin, tag)
        [fwhm_x_sted, fwhm_y_sted, fwhm_unfilt_sted, amp_sted, bg_sted, ROIgoodnessinfo_sted, halfROIused]...
            = fwhm_fit(cdata_sted, map_sted, ROIspanmax, pixelsize, rois_sted, mingoodness, FWHMmax, FWHMmin, tag, nroi2compare);
        numbeads_sted(ii) = length(fwhm_x_sted);
        fwhm_x_all_sted(1:numbeads_sted(ii),ii) = fwhm_x_sted;
        fwhm_y_all_sted(1:numbeads_sted(ii),ii) = fwhm_y_sted;
        fwhm_xyav_all_sted(1:numbeads_sted(ii),ii) = (fwhm_x_sted +fwhm_y_sted)/2;
        median_sted_all(ii,1) = median(fwhm_x_sted);
        median_sted_all(ii,2) = median(fwhm_y_sted);
        median_sted_all(ii,3) = median((fwhm_x_sted + fwhm_y_sted)/2); %take mean of x and y then find median of them
        mean_sted_all(ii,1) = mean(fwhm_x_sted);
        mean_sted_all(ii,2) = mean(fwhm_y_sted,1);
        mean_sted_all(ii,3) = mean((fwhm_x_sted + fwhm_y_sted)/2);
        std_sted_all(ii,1) = std(fwhm_x_sted);
        std_sted_all(ii,2) = std(fwhm_y_sted);
        std_sted_all(ii,3) = std((fwhm_x_sted + fwhm_y_sted)/2);
        amp_sted_all(1:numbeads_sted(ii),ii) = amp_sted;
        bg_sted_all(1:numbeads_sted(ii),ii) = bg_sted;
        ROIgoodness_all_sted(:, :,ii) = ROIgoodnessinfo_sted;
        fwhm_unfilt_sted_all(1:length(fwhm_unfilt_sted),:,ii) = fwhm_unfilt_sted(1:length(fwhm_unfilt_sted),:);
        STEDhalfROIused(1,ii) = halfROIused;
    elseif rois_sted == 0
        choice = 0;
    end

    if rois_dymin13 ~= 0
        %[fwhm_x, fwhm_y,fwhm_unfilt, amp, bg, ROIgoodnessinfo] =fwhm_fit(cdata, ...
        % ~, ROIspan, pixelsize, rois, mingoodness, FWHMmax, FWHMmin, tag)
        [fwhm_x_dymin13, fwhm_y_dymin13, fwhm_unfilt_dymin13, amp_dymin13, bg_dymin13, ROIgoodnessinfo_dymin13, halfROIused]...
            = fwhm_fit(cdata_dymin13, map_dymin13, ROIspanmax, pixelsize, rois_dymin13, mingoodness, FWHMmax, FWHMmin, tag, nroi2compare);
        numbeads_dymin13(ii) = length(fwhm_x_dymin13);
        fwhm_x_all_dymin13(1:numbeads_dymin13(ii),ii) = fwhm_x_dymin13;
        fwhm_y_all_dymin13(1:numbeads_dymin13(ii),ii) = fwhm_y_dymin13;
        fwhm_xyav_all_dymin13(1:numbeads_dymin13(ii),ii) = (fwhm_x_dymin13 + fwhm_y_dymin13)/2;
        median_dymin13_all(ii,1) = median(fwhm_x_dymin13);
        median_dymin13_all(ii,2) = median(fwhm_y_dymin13);
        median_dymin13_all(ii,3) = median((fwhm_x_dymin13 + fwhm_y_dymin13)/2); %take mean of x and y then find median of them
        mean_dymin13_all(ii,1) = mean(fwhm_x_dymin13);
        mean_dymin13_all(ii,2) = mean(fwhm_y_dymin13);
        mean_dymin13_all(ii,3) = mean((fwhm_x_dymin13 + fwhm_y_dymin13)/2);
        std_dymin13_all(ii,1) = std(fwhm_x_dymin13);
        std_dymin13_all(ii,2) = std(fwhm_y_dymin13);
        std_dymin13_all(ii,3) = std((fwhm_x_dymin13 + fwhm_y_dymin13)/2);
        amp_dymin13_all(1:numbeads_dymin13(ii),ii) = amp_dymin13;
        bg_dymin13_all(1:numbeads_dymin13(ii),ii) = bg_dymin13;
        ROIgoodness_all_dymin13(:, :,ii) = ROIgoodnessinfo_dymin13;
        fwhm_unfilt_dymin13_all(1:length(fwhm_unfilt_dymin13),:,ii) = fwhm_unfilt_dymin13(1:length(fwhm_unfilt_dymin13),:);
        DyMIN13halfROIused(1,ii) = halfROIused;
    elseif rois_dymin13 == 0
        choice = 0;
    end
    if rois_dymin14 ~= 0
        %[fwhm_x, fwhm_y,fwhm_unfilt, amp, bg, ROIgoodnessinfo] =fwhm_fit(cdata, ...
        % ~, ROIspan, pixelsize, rois, mingoodness, FWHMmax, FWHMmin, tag)
        [fwhm_x_dymin14, fwhm_y_dymin14, fwhm_unfilt_dymin14, amp_dymin14, bg_dymin14, ROIgoodnessinfo_dymin14, halfROIused]...
            = fwhm_fit(cdata_dymin14, map_dymin14, ROIspanmax, pixelsize, rois_dymin14, mingoodness, FWHMmax, FWHMmin, tag, nroi2compare);
        numbeads_dymin14(ii) = length(fwhm_x_dymin14);
        fwhm_x_all_dymin14(1:numbeads_dymin14(ii),ii) = fwhm_x_dymin14;
        fwhm_y_all_dymin14(1:numbeads_dymin14(ii),ii) = fwhm_y_dymin14;
        fwhm_xyav_all_dymin14(1:numbeads_dymin14(ii),ii) = (fwhm_x_dymin14 + fwhm_y_dymin14)/2;
        median_dymin14_all(ii,1) = median(fwhm_x_dymin14);
        median_dymin14_all(ii,2) = median(fwhm_y_dymin14);
        median_dymin14_all(ii,3) = median((fwhm_x_dymin14 + fwhm_y_dymin14)/2); %take mean of x and y then find median of them
        mean_dymin14_all(ii,1) = mean(fwhm_x_dymin14);
        mean_dymin14_all(ii,2) = mean(fwhm_y_dymin14);
        mean_dymin14_all(ii,3) = mean((fwhm_x_dymin14 + fwhm_y_dymin14)/2);
        std_dymin14_all(ii,1) = std(fwhm_x_dymin14);
        std_dymin14_all(ii,2) = std(fwhm_y_dymin14);
        std_dymin14_all(ii,3) = std((fwhm_x_dymin14 + fwhm_y_dymin14)/2);
        amp_dymin14_all(1:numbeads_dymin14(ii),ii) = amp_dymin14;
        bg_dymin14_all(1:numbeads_dymin14(ii),ii) = bg_dymin14;
        ROIgoodness_all_dymin14(:, :,ii) = ROIgoodnessinfo_dymin14;
        fwhm_unfilt_dymin14_all(1:length(fwhm_unfilt_dymin14),:,ii) = fwhm_unfilt_dymin14(1:length(fwhm_unfilt_dymin14),:);
        DyMIN14halfROIused(1,ii) = halfROIused;
    elseif rois_dymin14 == 0
        choice = 0;
    end
    if rois_dymin15 ~= 0
        %[fwhm_x, fwhm_y,fwhm_unfilt, amp, bg, ROIgoodnessinfo] =fwhm_fit(cdata, ...
        % ~, ROIspan, pixelsize, rois, mingoodness, FWHMmax, FWHMmin, tag)
        [fwhm_x_dymin15, fwhm_y_dymin15, fwhm_unfilt_dymin15, amp_dymin15, bg_dymin15, ROIgoodnessinfo_dymin15, halfROIused]...
            = fwhm_fit(cdata_dymin15, map_dymin15, ROIspanmax, pixelsize, rois_dymin15, mingoodness, FWHMmax, FWHMmin, tag, nroi2compare);
        numbeads_dymin15(ii) = length(fwhm_x_dymin15);
        fwhm_x_all_dymin15(1:numbeads_dymin15(ii),ii) = fwhm_x_dymin15;
        fwhm_y_all_dymin15(1:numbeads_dymin15(ii),ii) = fwhm_y_dymin15;
        fwhm_xyav_all_dymin15(1:numbeads_dymin15(ii),ii) = (fwhm_x_dymin15 + fwhm_y_dymin15)/2;
        median_dymin15_all(ii,1) = median(fwhm_x_dymin15);
        median_dymin15_all(ii,2) = median(fwhm_y_dymin15);
        median_dymin15_all(ii,3) = median((fwhm_x_dymin15 + fwhm_y_dymin15)/2); %take mean of x and y then find median of them
        mean_dymin15_all(ii,1) = mean(fwhm_x_dymin15);
        mean_dymin15_all(ii,2) = mean(fwhm_y_dymin15);
        mean_dymin15_all(ii,3) = mean((fwhm_x_dymin15 + fwhm_y_dymin15)/2);
        std_dymin15_all(ii,1) = std(fwhm_x_dymin15);
        std_dymin15_all(ii,2) = std(fwhm_y_dymin15);
        std_dymin15_all(ii,3) = std((fwhm_x_dymin15 + fwhm_y_dymin15)/2);
        amp_dymin15_all(1:numbeads_dymin15(ii),ii) = amp_dymin15;
        bg_dymin15_all(1:numbeads_dymin15(ii),ii) = bg_dymin15;
        ROIgoodness_all_dymin15(:, :,ii) = ROIgoodnessinfo_dymin15;
        fwhm_unfilt_dymin15_all(1:length(fwhm_unfilt_dymin15),:,ii) = fwhm_unfilt_dymin15(1:length(fwhm_unfilt_dymin15),:);
        DyMIN15halfROIused(1,ii) = halfROIused;
    elseif rois_dymin15 == 0
        choice = 0;
    end
    if rois_dymin16 ~= 0
        %[fwhm_x, fwhm_y,fwhm_unfilt, amp, bg, ROIgoodnessinfo] =fwhm_fit(cdata, ...
        % ~, ROIspan, pixelsize, rois, mingoodness, FWHMmax, FWHMmin, tag)
        [fwhm_x_dymin16, fwhm_y_dymin16, fwhm_unfilt_dymin16, amp_dymin16, bg_dymin16, ROIgoodnessinfo_dymin16, halfROIused]...
            = fwhm_fit(cdata_dymin16, map_dymin16, ROIspanmax, pixelsize, rois_dymin16, mingoodness, FWHMmax, FWHMmin, tag, nroi2compare);
        numbeads_dymin16(ii) = length(fwhm_x_dymin16);
        fwhm_x_all_dymin16(1:numbeads_dymin16(ii),ii) = fwhm_x_dymin16;
        fwhm_y_all_dymin16(1:numbeads_dymin16(ii),ii) = fwhm_y_dymin16;
        fwhm_xyav_all_dymin16(1:numbeads_dymin16(ii),ii) = (fwhm_x_dymin16 + fwhm_y_dymin16)/2;
        median_dymin16_all(ii,1) = median(fwhm_x_dymin16);
        median_dymin16_all(ii,2) = median(fwhm_y_dymin16);
        median_dymin16_all(ii,3) = median((fwhm_x_dymin16 + fwhm_y_dymin16)/2); %take mean of x and y then find median of them
        mean_dymin16_all(ii,1) = mean(fwhm_x_dymin16);
        mean_dymin16_all(ii,2) = mean(fwhm_y_dymin16);
        mean_dymin16_all(ii,3) = mean((fwhm_x_dymin16 + fwhm_y_dymin16)/2);
        std_dymin16_all(ii,1) = std(fwhm_x_dymin16);
        std_dymin16_all(ii,2) = std(fwhm_y_dymin16);
        std_dymin16_all(ii,3) = std((fwhm_x_dymin16 + fwhm_y_dymin16)/2);
        amp_dymin16_all(1:numbeads_dymin16(ii),ii) = amp_dymin16;
        bg_dymin16_all(1:numbeads_dymin16(ii),ii) = bg_dymin16;
        ROIgoodness_all_dymin16(:, :,ii) = ROIgoodnessinfo_dymin16;
        fwhm_unfilt_dymin16_all(1:length(fwhm_unfilt_dymin16),:,ii) = fwhm_unfilt_dymin16(1:length(fwhm_unfilt_dymin16),:);
        DyMIN16halfROIused(1,ii) = halfROIused;
    elseif rois_dymin16 == 0
        choice = 0;
    end
end
%% Replace zeros with NaNs so they won't be used in calculations
if choice ~= 0
    amp_sted_all(amp_sted_all==0)=NaN;
    bg_sted_all(bg_sted_all==0)=NaN;
    mean_sted_all(mean_sted_all==0)=NaN;
    median_sted_all(median_sted_all==0)=NaN;
    std_sted_all(std_sted_all==0)=NaN;

    fwhm_x_all_dymin13(fwhm_x_all_dymin13==0)=NaN;
    fwhm_y_all_dymin13(fwhm_y_all_dymin13==0)=NaN;
    fwhm_xyav_all_dymin13(fwhm_xyav_all_dymin13==0)=NaN;
    amp_dymin13_all(amp_dymin13_all==0)=NaN;
    bg_dymin13_all(bg_dymin13_all==0)=NaN;
    mean_dymin13_all(mean_dymin13_all==0)=NaN;
    median_dymin13_all(median_dymin13_all==0)=NaN;
    std_dymin13_all(std_dymin13_all==0)=NaN;

    fwhm_x_all_dymin14(fwhm_x_all_dymin14==0)=NaN;
    fwhm_y_all_dymin14(fwhm_y_all_dymin14==0)=NaN;
    fwhm_xyav_all_dymin14(fwhm_xyav_all_dymin14==0)=NaN;
    amp_dymin14_all(amp_dymin14_all==0)=NaN;
    bg_dymin14_all(bg_dymin14_all==0)=NaN;
    mean_dymin14_all(mean_dymin14_all==0)=NaN;
    median_dymin14_all(median_dymin14_all==0)=NaN;
    std_dymin14_all(std_dymin14_all==0)=NaN;

    fwhm_x_all_dymin15(fwhm_x_all_dymin15==0)=NaN;
    fwhm_y_all_dymin15(fwhm_y_all_dymin15==0)=NaN;
    fwhm_xyav_all_dymin15(fwhm_xyav_all_dymin15==0)=NaN;
    amp_dymin15_all(amp_dymin15_all==0)=NaN;
    bg_dymin15_all(bg_dymin15_all==0)=NaN;
    mean_dymin15_all(mean_dymin15_all==0)=NaN;
    median_dymin15_all(median_dymin15_all==0)=NaN;
    std_dymin15_all(std_dymin15_all==0)=NaN;

    fwhm_x_all_dymin16(fwhm_x_all_dymin16==0)=NaN;
    fwhm_y_all_dymin16(fwhm_y_all_dymin16==0)=NaN;
    fwhm_xyav_all_dymin16(fwhm_xyav_all_dymin16==0)=NaN;
    amp_dymin16_all(amp_dymin16_all==0)=NaN;
    bg_dymin16_all(bg_dymin16_all==0)=NaN;
    mean_dymin16_all(mean_dymin16_all==0)=NaN;
    median_dymin16_all(median_dymin16_all==0)=NaN;
    std_dymin16_all(std_dymin16_all==0)=NaN;

%     FWHM_all(:,1) = fwhm_x_all_sted*pixelsize;
%     FWHM_all(:,2) = fwhm_y_all_sted*pixelsize;
% 
%     stdevx = std(fwhm_x_all, 0, 2, 'omitnan')*pixelsize;
%     stdevy = std(fwhm_y_all, 0, 2, 'omitnan')*pixelsize;

% Plotting s
%     corrfactorx = 1.0103;
%     corrfactory = .9865;
%     figure
%     boxplot([FWHM_all(:,1)*corrfactorx, FWHM_all(:,2)*corrfactory], 'Labels', {'FWHM_x', 'FWHM_y'})
%     title(['# of spots: ', num2str(numbeads),' size limits = ',...
%         num2str(sigma_to_fwhm*pixelsize*FWHMmin), ' to ', num2str(sigma_to_fwhm*pixelsize*FWHMmax)])
%% 

    %edges2 = linspace(0, 400, 40);
%     figure(5)
%     for ii=1:num2compare
%         histogram((fwhm_x_all_sted(:,ii) + fwhm_y_all_sted(:,ii))/2*1000, 15)
%         hold on;
% %         histogram((fwhm_x_all(2,:) + fwhm_y_all(2,:))/2*1000, edges2)
%         % title('Histogram of FWHM fits (no filtering)')
% %         legend('2P','STED')
% 
%         title('Conv STED')
% %         xlim([40 70])
%         set(gca, 'fontsize', 18)
%     end
%    figure(6)
%     for ii=1:num2compare
%         histogram((fwhm_x_all_dymin13(:,ii) + fwhm_y_all_dymin13(:,ii))/2*1000, 15)
%         hold on;
% %         histogram((fwhm_x_all(2,:) + fwhm_y_all(2,:))/2*1000, edges2)
%         % title('Histogram of FWHM fits (no filtering)')
% %         legend('2P','STED')
% 
%         title('DyMIN')
% %         xlim([40 70])
%         set(gca, 'fontsize', 18)
%     end
%% Plot average Amplitude vs image # to get a nice photobleaching plot
figure(7)
% Amplitude over BG plot looked crazy.  Try plotting mean amplitude vs
% image #.
% amp_div_bg_sted_all=rmmissing(amp_sted_all(:,:))./rmmissing(bg_sted_all(:,:));
% amp_div_bg_sted = mean(amp_div_bg_sted_all,1);
% amp_div_bg_sted_std = std(amp_div_bg_sted_all,1);
% errorbar(1:num2compare,amp_div_bg_sted,amp_div_bg_sted_std);
% xlabel("Image Number")
% ylabel("Mean Amplitude/Background of beads in Image")
% amp_sted_mean = mean(rmmissing(amp_sted_all,1));
% amp_sted_std =std(rmmissing(amp_sted_all,1));
% amp_dymin16_mean = mean(rmmissing(amp_dymin16_all,1));
% amp_dymin16_std =std(rmmissing(amp_dymin16_all,1));
% amp_dymin15_mean = mean(rmmissing(amp_dymin15_all,1));
% amp_dymin15_std =std(rmmissing(amp_dymin15_all,1));
% amp_dymin14_mean = mean(rmmissing(amp_dymin14_all,1));
% amp_dymin14_std =std(rmmissing(amp_dymin14_all,1));
% amp_dymin13_mean = mean(rmmissing(amp_dymin13_all,1));
% amp_dymin13_std =std(rmmissing(amp_dymin13_all,1));

amp_sted_mean = mean(amp_sted_all,'omitnan');
amp_sted_median = median(amp_sted_all,'omitnan');
amp_sted_std = std(amp_sted_all,'omitnan');

amp_dymin13_mean = mean(amp_dymin13_all,'omitnan');
amp_dymin13_median = median(amp_dymin13_all,'omitnan');
amp_dymin13_std = std(amp_dymin13_all,'omitnan');

amp_dymin14_mean = mean(amp_dymin14_all,'omitnan');
amp_dymin14_median = median(amp_dymin14_all,'omitnan');
amp_dymin14_std = std(amp_dymin14_all,'omitnan');

amp_dymin15_mean = mean(amp_dymin15_all,'omitnan');
amp_dymin15_median = median(amp_dymin15_all,'omitnan');
amp_dymin15_std = std(amp_dymin15_all,'omitnan');

amp_dymin16_mean = mean(amp_dymin16_all,'omitnan');
amp_dymin16_median = median(amp_dymin16_all,'omitnan');
amp_dymin16_std = std(amp_dymin16_all,'omitnan');
%errorbar(1:num2compare,amp_sted_mean,amp_sted_std)
plot(1:num2compare,amp_sted_mean)
hold on
%errorbar(1:num2compare,amp_dymin_mean,amp_dymin_std)
plot(1:num2compare,amp_dymin16_mean)
plot(1:num2compare,amp_dymin15_mean)
plot(1:num2compare,amp_dymin14_mean)
plot(1:num2compare,amp_dymin13_mean)
set(gca, 'fontsize', 18)
xlabel("Image Number")
ylabel("Mean Amplitude of beads in Image")
legend('Conv STED', 'DyMIN1','DyMIN2','DyMIN3','DyMIN4')
title('Photobleaching over time')
%% Plot mean with std error bars for FWHM for conv STED, x, y and avg
figure(8)
plot([1:50],mean_sted_all(:,3)*1000)
hold on
plot([1:50],mean_sted_all(:,1)*1000)
plot([1:50],mean_sted_all(:,2)*1000)
set(gca, 'fontsize', 18)
legend("avg of x and y","x","y")
xlabel("Image Number")
ylabel("Mean FWHM of Lorentzian fits of beads, ConvSTED (nm)")
title("Conv STED")

%% Plot mean with std error bars for FWHM, x, y and avg
figure(9)
subplot(2,2,1)
plot([1:50],mean_dymin16_all(:,3)*1000)
hold on
plot([1:50],mean_dymin16_all(:,1)*1000)
plot([1:50],mean_dymin16_all(:,2)*1000)
set(gca, 'fontsize', 18)
legend("avg of x and y","x","y")
xlabel("Image Number")
ylabel("Mean FWHM of Lorentzian fits of beads(nm)")
title("DyMIN1")
subplot(2,2,2)
plot([1:50],mean_dymin15_all(:,3)*1000)
hold on
plot([1:50],mean_dymin15_all(:,1)*1000)
plot([1:50],mean_dymin15_all(:,2)*1000)
set(gca, 'fontsize', 18)
legend("avg of x and y","x","y")
xlabel("Image Number")
ylabel("Mean FWHM of Lorentzian fits of beads(nm)")
title("DyMIN2")
subplot(2,2,3)
plot([1:50],mean_dymin14_all(:,3)*1000)
hold on
plot([1:50],mean_dymin14_all(:,1)*1000)
plot([1:50],mean_dymin14_all(:,2)*1000)
set(gca, 'fontsize', 18)
legend("avg of x and y","x","y")
xlabel("Image Number")
ylabel("Mean FWHM of Lorentzian fits of beads(nm)")
title("DyMIN3")
subplot(2,2,4)
plot([1:50],mean_dymin13_all(:,3)*1000)
hold on
plot([1:50],mean_dymin13_all(:,1)*1000)
plot([1:50],mean_dymin13_all(:,2)*1000)
set(gca, 'fontsize', 18)
legend("avg of x and y","x","y")
xlabel("Image Number")
ylabel("Mean FWHM of Lorentzian fits of beads(nm)")
title("DyMIN4")
%% Plot mean with std error bars for FWHM with just conv STED and DyMIN
figure(10)
plot([1:50],mean_sted_all(:,3)*1000)
hold on
plot([1:50],mean_dymin16_all(:,3)*1000)
plot([1:50],mean_dymin15_all(:,3)*1000)
plot([1:50],mean_dymin14_all(:,3)*1000)
plot([1:50],mean_dymin13_all(:,3)*1000)
set(gca, 'fontsize', 18)
legend('Conv. STED', 'DyMIN1','DyMIN2','DyMIN3','DyMIN4')
title('Conv STED vs DyMIN mean bead fit FWHM, average of x and y')
xlabel("Image Number")
ylabel("Mean FWHM of Lorentzian fits of beads (nm)")

%% Plot median with std error bars for FWHM with just conv STED and DyMIN
figure(11)
plot([1:50],median_sted_all(:,3)*1000)
hold on
plot([1:50],median_dymin16_all(:,3)*1000)
plot([1:50],median_dymin15_all(:,3)*1000)
plot([1:50],median_dymin14_all(:,3)*1000)
plot([1:50],median_dymin13_all(:,3)*1000)
set(gca, 'fontsize', 18)
legend('Conv. STED', 'DyMIN1','DyMIN2','DyMIN3','DyMIN4')
title('Conv STED vs DyMIN median bead fit FWHM, average of x and y')
xlabel("Image Number")
ylabel("Median FWHM of Lorentzian fits of beads (nm)")
%% Make boxplots of bead fit FWHM vs image #.  Compare x, y and avg for Conv STED
figure(12)
boxplot(fwhm_x_all_sted, 'Colors','b')
hold on
boxplot(fwhm_y_all_sted,'Colors', 'r')
boxplot(fwhm_xyav_all_sted,'Colors','g')
title('Conv STED')


%% Make boxplots of bead fit FWHM vs image #. Compare x, y and avg for DyMIN
figure(13)
boxplot(fwhm_x_all_dymin13, 'Colors','b')
hold on
boxplot(fwhm_y_all_dymin13,'Colors', 'r')
boxplot(fwhm_xyav_all_dymin13,'Colors','g')
title('DyMIN4')
%% Make boxplots of bead fit FWHM vs image #. Compare Conv STED and DyMIN
figure(14)
boxplot(fwhm_xyav_all_sted, 'Colors','b','Widths',2)
hold on
boxplot(fwhm_xyav_all_dymin13,'Colors', 'r','Widths',2)
set(gca, 'fontsize', 18)
title('Conv STED blue, DyMIN red')
%         avg_unfilt_nosted = (fwhm_unfilt_nosted(:,1)+fwhm_unfilt_nosted(:,2))/2;
%         avg_unfilt_sted = (fwhm_unfilt_convSTED(:,1)+fwhm_unfilt_convSTED(:,2))/2;
%      for ii=1:num2compare
%         figure
%         histogram(avg_unfilt_nosted(, edges2)
%         hold on;
%         histogram(avg_unfilt_sted, edges2)
%         % title('Histogram of FWHM fits (no filtering)')
%         legend('2P','STED')
%         ylabel('Number of beads')
%         xlabel('FWHM of fit (nm)')
%         title('Raw distribution')
%         xlim([0 inf])
%         set(gca, 'fontsize', 18)
%      end 
        %%
%         median_sted_unfilt = median(avg_unfilt_sted)
  %      median_sted_filt = (fwhm_x_median(2) + fwhm_y_median(2))/2
%         median_nosted_unfilt = median(avg_unfilt_nosted)
 %       median_nosted_filt = (fwhm_x_median(1) + fwhm_y_median(1))/2
        
        %% for res vs power measurements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         powers = .51*[0 13 22.4 34 44.7 55.5 66 79];
%         x0 = [1, 2];
%         theory = lsqcurvefit(@theo_res_fit, x0, powers, fwhm_x_avg);
%         fit = theo_res_fit(theory, powers);
%         figure
%         ax = gca;
%         boxplot(fwhm_x_all'*pixelsize, 'colors', 'b', 'PlotStyle', 'traditional', 'OutlierSize', .01)
%         hold on;
%         boxplot(fwhm_y_all'*pixelsize, 'colors', 'k', 'PlotStyle', 'traditional',  'OutlierSize', .01)
%         hold on;
%         plot(fit, 'linewidth', 4)
%         ax.XTick = 1:num2compare;
%         ax.XTickLabel = (num2str(powers', '%.2f'));
%         ax.FontSize = 12;
%         xlabel('STED power at objective (mW)')
%         ylabel('FHWM (\mum)')
%         ylim([.08 .32])
%         title(['# of spots: ', num2str(numbeads)])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Making plots for DyMIN v1 paper
xx = [1:0.01:num2compare]';
figure(15)
% subplot(1,2,1)
% plot(1:num2compare,amp_sted_median/mean(amp_sted_median(1:3)),'*r')
% hold on
% 
% plot(1:num2compare,amp_dymin16_median/mean(amp_dymin16_median(1:3)),'+g')
% 
% plot(1:num2compare,amp_dymin15_median/mean(amp_dymin15_median(1:3)),'ob')
% 
% plot(1:num2compare,amp_dymin14_median/mean(amp_dymin14_median(1:3)),'.c')
% 
% plot(1:num2compare,amp_dymin13_median/mean(amp_dymin13_median(1:3)),'xm')
% s0 = spline(1:num2compare,amp_sted_median/mean(amp_sted_median(1:3)),xx);
% s1 = spline(1:num2compare,amp_dymin16_median/mean(amp_dymin16_median(1:3)),xx);
% s2 = spline(1:num2compare,amp_dymin15_median/mean(amp_dymin15_median(1:3)),xx);
% s3 = spline(1:num2compare,amp_dymin14_median/mean(amp_dymin14_median(1:3)),xx);
% s4 = spline(1:num2compare,amp_dymin13_median/mean(amp_dymin13_median(1:3)),xx);
% 
% plot(xx,s0,'r');
% plot(xx,s1,'g');
% plot(xx,s2,'b');
% plot(xx,s3,'c');
% plot(xx,s4,'m');
% %set(gca, 'fontsize', 18)
% xlabel("Image Number")
% ylabel("Median Amplitude scaled at t=0")
% legend('Conv STED', 'DyMIN1','DyMIN2','DyMIN3','DyMIN4')
% plot(1:num2compare,amp_sted_mean)
% hold on
% %errorbar(1:num2compare,amp_dymin_mean,amp_dymin_std)
% plot(1:num2compare,amp_dymin16_mean)
% plot(1:num2compare,amp_dymin15_mean)
% plot(1:num2compare,amp_dymin14_mean)
% plot(1:num2compare,amp_dymin13_mean)
% set(gca, 'fontsize', 18)
% xlabel("Image Number")
% ylabel("Mean Amplitude of beads in Image")
% legend('Conv STED', 'DyMIN1','DyMIN2','DyMIN3','DyMIN4')
% title('Photobleaching over time')
subplot(1,2,1)

plot(1:num2compare,amp_dymin16_mean/mean(amp_dymin16_mean(1:3)),'+g')
hold on
plot(1:num2compare,amp_dymin15_mean/mean(amp_dymin15_mean(1:3)),'ob')

plot(1:num2compare,amp_dymin14_mean/mean(amp_dymin14_mean(1:3)),'.c')

plot(1:num2compare,amp_dymin13_mean/mean(amp_dymin13_mean(1:3)),'xm')
plot(1:num2compare,amp_sted_mean/mean(amp_sted_mean(1:3)),'*r')

% s0 = spline(1:num2compare,amp_sted_mean/mean(amp_sted_mean(1:3)),xx);
% s1 = spline(1:num2compare,amp_dymin16_mean/mean(amp_dymin16_mean(1:3)),xx);
% s2 = spline(1:num2compare,amp_dymin15_mean/mean(amp_dymin15_mean(1:3)),xx);
% s3 = spline(1:num2compare,amp_dymin14_mean/mean(amp_dymin14_mean(1:3)),xx);
% s4 = spline(1:num2compare,amp_dymin13_mean/mean(amp_dymin13_mean(1:3)),xx);
% plot(xx,s0,'r');
% plot(xx,s1,'g');
% plot(xx,s2,'b');
% plot(xx,s3,'c');
% plot(xx,s4,'m');
%set(gca, 'fontsize', 18)
xlabel("Image Number")
axis([1 num2compare 0 1.05])
ylabel("Mean Amplitude scaled at t=0")
legend( 'DyMIN1','DyMIN2','DyMIN3','DyMIN4','STED')
subplot(1,2,2)
plot(1:num2compare,numbeads_dymin16,'+g')
hold on
plot(1:num2compare,numbeads_dymin15,'ob')
plot(1:num2compare,numbeads_dymin14,'.c')
plot(1:num2compare,numbeads_dymin13,'xm')
plot(1:num2compare,numbeads_sted,'*r')
% s0 = spline(1:num2compare,numbeads_sted,xx);
% s1 = spline(1:num2compare,numbeads_dymin16,xx);
% s2 = spline(1:num2compare,numbeads_dymin15,xx);
% s3 = spline(1:num2compare,numbeads_dymin14,xx);
% s4 = spline(1:num2compare,numbeads_dymin13,xx);
% plot(xx,s0,'r');
% plot(xx,s1,'g');
% plot(xx,s2,'b');
% plot(xx,s3,'c');
% plot(xx,s4,'m');
xlabel("Image Number")
axis([1 num2compare 0 215]);
ylabel("Number of Beads")




end

%% Make bar graph of after/before confocal ROI counts ratio
figure(16)
con_amp_ratio = [0.446 0.387 0.454 0.357 0.310];
X = categorical({'DyMIN1','DyMIN2','DyMIN3','DyMIN4','STED'});
X = reordercats(X,{'DyMIN1','DyMIN2','DyMIN3','DyMIN4','STED'});
Y = con_amp_ratio;
b = bar(X,Y);
b.FaceColor = 'flat';
b.CData(5,:) = [0.5 0 0.5];
ylabel("Counts in Confocal ROI after/before")
%% Plot mean bead amplitudevs image number w 4th degree poly fit
fig17=figure(17);


plot(1:num2compare,amp_dymin16_mean/mean(amp_dymin16_mean(1:3)),'+g')
hold on
plot(1:num2compare,amp_dymin15_mean/mean(amp_dymin15_mean(1:3)),'ob')
plot(1:num2compare,amp_dymin14_mean/mean(amp_dymin14_mean(1:3)),'.c')
plot(1:num2compare,amp_dymin13_mean/mean(amp_dymin13_mean(1:3)),'xm')
plot(1:num2compare,amp_sted_mean/mean(amp_sted_mean(1:3)),'*r')


%Do poly fits to data for lines, order of 4 was what we agreed matched the data best
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
fig18=figure(18);
plot(1:num2compare,numbeads_dymin16/mean(numbeads_dymin16(1:3)),'+g')
hold on
plot(1:num2compare,numbeads_dymin15/mean(numbeads_dymin15(1:3)),'ob')
plot(1:num2compare,numbeads_dymin14/mean(numbeads_dymin14(1:3)),'.c')
plot(1:num2compare,numbeads_dymin13/mean(numbeads_dymin13(1:3)),'xm')
plot(1:num2compare,numbeads_sted/mean(numbeads_sted(1:3)),'*r')


%Do poly fits to data for lines
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
axis([1 num2compare 50 120])

% fig18.Units               = 'centimeters';
% fig18.Position(3)         = 8;
% fig18.Position(4)         = 6;
% set(fig18.Children, ...
%     'FontName',     'Times', ...
%     'FontSize',     9);
% set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
% fig18.PaperPositionMode   = 'auto';
% print('Q:\OneDrive - The University of Colorado Denver\From Dropbox\STED\STED Data\230822\45nmYG_DyMINvsSTED\Numbeads_Poly', '-dpng', '-r600')


%% Maybe just show the data???
figure(19)
subplot(1,2,1)

plot(1:num2compare,amp_dymin16_mean/mean(amp_dymin16_mean(1:3)),'+g')
hold on
plot(1:num2compare,amp_dymin15_mean/mean(amp_dymin15_mean(1:3)),'ob')

plot(1:num2compare,amp_dymin14_mean/mean(amp_dymin14_mean(1:3)),'.c')

plot(1:num2compare,amp_dymin13_mean/mean(amp_dymin13_mean(1:3)),'xm')
plot(1:num2compare,amp_sted_mean/mean(amp_sted_mean(1:3)),'*r')
xlabel("Image Number")
axis([1 num2compare 0 1.05])
ylabel("Mean Amplitude scaled at t=0")
legend( 'DyMIN1','DyMIN2','DyMIN3','DyMIN4','STED')
subplot(1,2,2)
plot(1:num2compare,numbeads_dymin16/mean(numbeads_dymin16(1:3)),'+g')
hold on
plot(1:num2compare,numbeads_dymin15/mean(numbeads_dymin15(1:3)),'ob')
plot(1:num2compare,numbeads_dymin14/mean(numbeads_dymin14(1:3)),'.c')
plot(1:num2compare,numbeads_dymin13/mean(numbeads_dymin13(1:3)),'xm')
plot(1:num2compare,numbeads_sted/mean(numbeads_sted(1:3)),'*r')
xlabel("Image Number")
axis([1 num2compare 0 1.2]);
ylabel("Number of Beads scaled at t=0")