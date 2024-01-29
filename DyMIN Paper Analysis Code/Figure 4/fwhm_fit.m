function [fwhm_x, fwhm_y, fwhm_x_median, fwhm_y_median, fwhm_unfilt, ROIgoodnessinfo] =...
    fwhm_fit(cdata, ~, ROIspan, pixelsize, rois, mingoodness, FWHMmax, FWHMmin, tag,numROIs)
%%  rois is a list of x and y values for the ROI centers from find_peaks
%% Extract FWHM
halfROImax = (ROIspan-1)/2;
%numROIs = 11;

% parameters are: [Amplitude, x0, sigmax, y0, sigmay, background]
x0 = [10, 0, .5, 0, .5, 5]; %Inital guess parameters
%x = [2,2.2,7,3.4,4.5,+0.02*2*pi]; %centroid parameters
sigma_2_fwhm = 2*sqrt(2*log(2));%convert gaussian 1/e to fwhm
maxCount = double(max(max(cdata)));

InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'

clims = [1 maxCount/2];
%list of indicies for data in the 'fits' variable that don't meet certain
%specifications, like goodness of fit or unrealistic size of fit
outlier = [];

%toggle the usage of the lorentzian fit
if 1 == strcmp('sted', tag)
%     lorFit = 0;%use a gaussian fit
    lorFit = 1;
else
    lorFit = 0;
end

%% loop through all ROI centers, define a 2d patch, and fit
testfits = zeros(6, length(rois), numROIs);
goodnessArray = zeros(length(rois), numROIs);
halfROI = halfROImax + 1 - (1:numROIs);

for jj = 1:numROIs
    for ii = 1:length(rois)
        [X,Y] = meshgrid(-halfROI(jj):halfROI(jj));
        xdata = zeros(size(X,1),size(Y,2),2);
        xdata(:,:,1) = X;
        xdata(:,:,2) = Y;

%         numptshr = 301;
%         [Xhr,Yhr] = meshgrid(linspace(-halfROI,halfROI, numptshr)); % generate high res grid for plot
%         xdatahr = zeros(numptshr, numptshr, 2);
%         xdatahr(:,:,1) = Xhr;
%         xdatahr(:,:,2) = Yhr;

        x = (rois(ii,1));
        y = (rois(ii,2));
        Z = double(cdata(y-halfROI(jj):y+halfROI(jj), x-halfROI(jj):x+halfROI(jj)));

        lb = [0, -halfROI(jj), 0, -halfROI(jj), 0, 0];
        ub = [maxCount, halfROI(jj), ROIspan, halfROI(jj), ROIspan, x0(1)];
        % returns a vector of 6 elements fit(1) is amplidude
        % fit(2) is x coord, fit(3) is x waist
        % fit(4) is y coord, fit(5) is y waist, and fit(6) is
        % background
        if lorFit == 1 
            [fits, resnorm, ~, ~] = lsqcurvefit(@D2LorentianFunction, x0, xdata, Z, lb, ub);
        else
            [fits, resnorm, ~, ~] = lsqcurvefit(@D2GaussFunction, x0, xdata, Z, lb, ub);
        end

        %goodnessoffit defined in the supplemental material for 2018 optica
        %paper from vicidomini group
        goodnessoffit = 1 - resnorm/sum(sum((Z-mean(mean(Z))).^2));
        testfits(:, ii, jj) = fits;
        goodnessArray(ii, jj) = goodnessoffit;
    end
        
%     %%%%%plot the raw data and the fit randomly output about 1/15 of plots
%     if randi(15) == 15
%         fac = pixel+1;
%         if lorFit == 1
%             fit = D2LorentianFunction(fits, xdata);
%         else
%             fit = D2GaussFunction(fits, xdata);
%         end
% 
%         figure
%         subplot(2, 1, 1)
%         errorbar(xdata(1,:,1), fit(fac+round(fits(4)), :), sqrt(fit(fac+round(fits(4)), :)))
%         hold on
%         plot(xdata(1,:,1), Z(fac+round(fits(4)),:))
%         hold off
%         legend('x fit', 'x cut')
%         title(['FWHM_x Norm Res = ',num2str(resnorm) ])
%         subplot(2, 1, 2)
%         plot(xdata(1,:,1), fit(:,fac+round(fits(2))), xdata(1,:,1), Z(:,fac+round(fits(2))))
%         legend('y fit', 'y cut')
%         title('FWHM_y')
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
goodnesssum = sum(goodnessArray);
[~, best] = max(goodnesssum);
    
%each column is an entry for the fit parameters
fits_all = testfits(:,:,best);
% exitflag_all(:,ii) = exitflag;
goodnessoffit_all = goodnessArray(:, best);
ROIgoodnessinfo = [2*halfROI+1; goodnesssum];

%% %%%%%%%%% convert pixel to microns and gaussian waist to FWHM
if lorFit == 1
    fwhm_unfilt(:,1) = fits_all(3,:)*pixelsize*1000;
    fwhm_unfilt(:,2) = fits_all(5,:)*pixelsize*1000;
else
    fwhm_unfilt(:,1) = sigma_2_fwhm*fits_all(3,:)*pixelsize*1000;
    fwhm_unfilt(:,2) = sigma_2_fwhm*fits_all(5,:)*pixelsize*1000;
end

if lorFit == 1
    fits_all(3,:) = pixelsize*fits_all(3,:);
    fits_all(5,:) = pixelsize*fits_all(5,:);
else
    fits_all(3,:) = sigma_2_fwhm*pixelsize*fits_all(3,:);
    fits_all(5,:) = sigma_2_fwhm*pixelsize*fits_all(5,:);
end

%plot goodness of fit
figure
histogram(goodnessoffit_all, 80)
title('Plot showing the goodness of fit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(fits_all)
    %the center of the fit along x
    fitvalx_all(ii) = rois(ii,1) + fits_all(2,ii);
    %the center of the fit along y
    fitvaly_all(ii) = rois(ii,2) + fits_all(4,ii);
%     plot(fitvalx_all(i),fitvaly_all(i),'g*')
end

%Save data without removing residuals
%following three variables store the ROI centers
fitvalx_all_orig = fitvalx_all;
fitvaly_all_orig = fitvaly_all;
fits_all_orig = fits_all;

%% remove fits that don't meet certain criteria
%add indices of bad fits from fits_all to a list of bad entires in outliers
%remove low goodness of fit
kk = 1;
for ii=1:length(fits_all)
    if goodnessoffit_all(ii) < mingoodness 
        outlier(kk) = ii;
        kk = kk + 1;
    end
end

%remove large width value along x or small widths 
for ii=1:length(fits_all)
    if fits_all(3,ii) > FWHMmax || fits_all(3,ii) < FWHMmin
        outlier(kk) = ii;
        kk = kk + 1;
    end
end

%remove large width value along y
for ii=1:length(fits_all)
    if fits_all(5,ii) > FWHMmax || fits_all(5,ii) < FWHMmin
        outlier(kk) = ii;
        kk = kk + 1;
    end
end

%%
outlier = unique(outlier);
fits_all(:,outlier)=[];

figure
imagesc(cdata,clims)
colormap(hot)
pbaspect([ 1 1 1])
hold on
%See which were removed
fitvalx_all(outlier)=[];
fitvaly_all(outlier)=[];
for ii=1:length(fitvalx_all)
    plot(fitvalx_all(ii),fitvaly_all(ii),'g*')
end


%% amplitude vs width plot
figure
% Compare amplitude and width - see no correlation
plot(fits_all(1,:),fits_all(3,:),'k*')
xlabel('counts')
ylabel('width (microns)')
title('Fit amplitude vs fit width (FWHM)')


%%
%%%%%%%%%%%% plotting the waist in x vs the waist in y %%%%%%%%%%%%%%%%%%
% edges = linspace(0, 10, 60);
% edges2 = linspace(0, 400, 40);
%%%%%% unfilt %%%%%%%%%%%%%%%%%%%%%%
% figure
% histogram(fwhm_unfilt(:,1), edges2)
% hold on;
% histogram(fwhm_unfilt(:,2), edges2)
% ylabel('Number of beads')
% xlabel('FWHM of Lorentzian fit (nm)')
% % xlim([0 inf])
% legend('FWHM_x', 'FWHM_y')
% title('Histogram of FWHM fits (unfiltered)')
% set(gca, 'fontsize', 16)

%%%%%% filt %%%%%%%%%%%%%%%%%%%%%%
% figure
% histogram(fits_all(3,:)*pixelsize*1000, edges2)
% hold on;
% histogram(fits_all(5,:)*pixelsize*1000, edges2)
% ylabel('Number of beads')
% xlabel('FWHM of Lorentzian fit (nm)')
% % xlim([0 inf])
% legend('FWHM_x', 'FWHM_y')
% title('Histogram of FWHM fits')
% set(gca, 'fontsize', 16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% % fwhm in x and y
% plot(fits_all(3,:)*pixelsize, fits_all(5,:)*pixelsize,'g*')
% xlabel('fwhm_x')
% ylabel('fwhm_y')


%% analyze
outlier_all = outlier;
outlier_all = unique(outlier_all);

% fits_all(:,outlier_all)=[];
% 
% fitvalx_all(outlier_all)=[];
% fitvaly_all(outlier_all)=[];

fwhm_x = fits_all(3,:);
fwhm_y = fits_all(5,:);

fwhm_x_median = median(fwhm_x);
fwhm_y_median = median(fwhm_y);
stdx = std(fwhm_x);
stdy = std(fwhm_y);
stdavg = (stdx + stdy)/2;

%% %%%%%%%%%%%% remove data outside of some  number of standard deviations
% numstds = 2;
% fwhm_x_upper = fwhm_x_median + numstds*stdavg;
% fwhm_x_lower = fwhm_x_median - numstds*stdavg;
% fwhm_y_upper = fwhm_y_median + numstds*stdavg;
% fwhm_y_lower = fwhm_y_median - numstds*stdavg;
% kk = 1;
% for ii=1:length(fwhm_x)
%     if fwhm_x(ii) > fwhm_x_upper || fwhm_x(ii) < fwhm_x_lower
%         outlier_std(kk) = ii;
%         kk = kk + 1;
%     end
% end
% 
% for ii=1:length(fwhm_y)
%     if fwhm_y(ii) > fwhm_y_upper || fwhm_y(ii) < fwhm_y_lower
%         outlier_std(kk) = ii;
%         kk = kk + 1;
%     end
% end
% 
% outlier_std = unique(outlier_std);
% fwhm_x(outlier_std)=[];
% fwhm_y(outlier_std)=[];


disp("The file run was fwhm_fit in resolution_estimator folder in MATLAB")