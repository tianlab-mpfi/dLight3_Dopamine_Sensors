% This script makes heatmaps of dff0 values for all ROIs and for all
% animals, scaled 

% NOTE: first run looming_roi_analysis_with_dff0 or
% tailshock_roi_analysis_with_dff0 to get the nROIs_dff0_ROIs.mat file

function looming_heatmap_all_ROIs(filePath)

close all
clear all


% filePath = 'E:\data_summaries\dLight\dff0_roi_looming_datasets\data used in figures';
filePath = 'E:\data_summaries\dLight\dff0_roi_tailshock_datasets';

turbo = 0;


% Create the result directory if it doesn't exist
result = 'heatmaps_scale_min_max_dff0';
if ~exist(fullfile(filePath,result))
    mkdir(fullfile(filePath,result))
end

% Load data for dff0
roi_dff0_all = [];
file = '20*';
files = dir(fullfile(filePath,file));
for ii = 1:numel(files)
    load(fullfile(filePath,files(ii).name));
    roi_dff0_all(:,:,ii) = dff0_rois;
end


min_dff0 = min(min(min(roi_dff0_all)));
max_dff0 = max(max(max(roi_dff0_all)));


% Create heatmap for each dataset
for ii = 1:numel(files)
    dff0_thisData = roi_dff0_all(:,:,ii);
    % Sort ROIs by the way they are numbered in the ROI grid, which is 
    % 1 2 3
    % 7 6 4
    % % 8 9 5
    dff0_thisData = dff0_thisData';
    dff0_sort = zeros(size(dff0_thisData,1),size(dff0_thisData,2));
    dff0_sort(1:3,:) = dff0_thisData(1:3,:);
    dff0_sort(4,:) = dff0_thisData(7,:);
    dff0_sort(5,:) = dff0_thisData(6,:);
    dff0_sort(6,:) = dff0_thisData(4,:);
    dff0_sort(7,:) = dff0_thisData(8,:);
    dff0_sort(8,:) = dff0_thisData(9,:);
    dff0_sort(9,:) = dff0_thisData(5,:);
    dff0_plot = dff0_sort;
    % Plot values
    heatmap = imagesc(dff0_plot);
    if turbo == 1
        colormap('turbo')
    end
    min_scale = round(min_dff0-0.05,1);
    max_scale = round(max_dff0+0.05,1);
    clim([min_scale max_scale]);
    % clim([-2.2 2.2]); %for looming
    % clim([-1.7 2.7]); %for tailshock
    colorbar;
    figname = num2str(files(ii).name);
    figname = erase(figname,'.mat'); %get rid of file extension
    imwrite(dff0_plot,fullfile(filePath,result,figname),'tif')    
    saveas(heatmap,fullfile(filePath,result,figname),'png')
    saveas(heatmap,fullfile(filePath,result,figname),'pdf')
    saveas(heatmap,fullfile(filePath,result,figname),'fig')
end






