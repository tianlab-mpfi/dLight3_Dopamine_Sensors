% This script does a paired comparison test of the average peak of trials from 
% early part of experiment vs peak of trials from late part of experiment 

close all
clear all

horiz_displace = 0; %indicate if you want the data of separate sensors plotted at different places along x axis.

filePath = 'E:\data_summaries\dLight\tailshock_early_vs_late_trials\separate_sensors\dLight3_only';
filePath_dLight3a6 = 'E:\data_summaries\dLight\tailshock_early_vs_late_trials\tailshock_early_vs_late_dLight3a6';
filePath_dLight3a8 = 'E:\data_summaries\dLight\tailshock_early_vs_late_trials\tailshock_early_vs_late_dLight3a8';
file_dLight3a6 = 'peak*';
file_dLight3a8 = 'peak*';

all_files_dLight3a6 = dir(fullfile(filePath_dLight3a6,file_dLight3a6));
all_files_dLight3a8 = dir(fullfile(filePath_dLight3a8,file_dLight3a8));


% Combine data for dLight3a6
early_peak_dLight3a6 = [];
late_peak_dLight3a6 = [];

for ii = 1:size(all_files_dLight3a6)
    load(fullfile(filePath_dLight3a6,all_files_dLight3a6(ii).name));
    early_peak_dLight3a6 = [early_peak_dLight3a6;early_trials_peak_allROI(:,1)];
    late_peak_dLight3a6 = [late_peak_dLight3a6;late_trials_peak_allROI(:,1)];
end

early_peak_dLight3a6_avg = mean(early_peak_dLight3a6);
late_peak_dLight3a6_avg = mean(late_peak_dLight3a6);
early_peak_dLight3a6_median = median(early_peak_dLight3a6);
late_peak_dLight3a6_median = median(late_peak_dLight3a6);


% Combine data for dLight3a8
early_peak_dLight3a8 = [];
late_peak_dLight3a8 = [];

for ii = 1:size(all_files_dLight3a8)
    load(fullfile(filePath_dLight3a8,all_files_dLight3a8(ii).name));
    early_peak_dLight3a8 = [early_peak_dLight3a8;early_trials_peak_allROI(:,1)];
    late_peak_dLight3a8 = [late_peak_dLight3a8;late_trials_peak_allROI(:,1)];
end

early_peak_dLight3a8_avg = mean(early_peak_dLight3a8);
late_peak_dLight3a8_avg = mean(late_peak_dLight3a8);
early_peak_dLight3a8_median = median(early_peak_dLight3a8);
late_peak_dLight3a8_median = median(late_peak_dLight3a8);

% Average values across all sensors
early_peak_all = [early_peak_dLight3a6; early_peak_dLight3a8];
late_peak_all = [late_peak_dLight3a6; late_peak_dLight3a8];
early_peak_all_avg = mean(early_peak_all);
late_peak_all_avg = mean(late_peak_all);

% Combine all datapoints to determine limits for y axis
all_data = [early_peak_all; late_peak_all];



% Statistical test: sign rank
disp('signrank')
% % Test hypothesis that the the distribution mean of early trials is smaller
% % than the mean of late trials (i.e. the response increases over trials)
p_response_incr_signrank_dLight3a6 = signrank(early_peak_dLight3a6,late_peak_dLight3a6,'tail','left');
p_response_incr_signrank_dLight3a8 = signrank(early_peak_dLight3a8,late_peak_dLight3a8,'tail','left');
% % Test hypothesis that the the distribution mean of early trials is larger
% % than the mean of late trials (i.e. the response decreases over trials)
p_response_decr_signrank_dLight3a6 = signrank(early_peak_dLight3a6,late_peak_dLight3a6,'tail','right')
p_response_decr_signrank_dLight3a8 = signrank(early_peak_dLight3a8,late_peak_dLight3a8,'tail','right')


% Statistical test: t test
% same mean)
disp('t test')
% Test hypothesis that the the distribution mean of early trials is less
% than the mean of late trials (i.e. the response DEcreases over trials, since this is a negative response)
[h,p_response_incr_ttest_dLight3a6] = ttest(early_peak_dLight3a6,late_peak_dLight3a6,'tail','left');
[h,p_response_incr_ttest_dLight3a8] = ttest(early_peak_dLight3a8,late_peak_dLight3a8,'tail','left');
% Test hypothesis that the the distribution mean of early trials is larger
% than the mean of late trials (i.e. the response INcreases over trials, since this is a negative response)
[h,p_response_decr_ttest_dLight3a6] = ttest(early_peak_dLight3a6,late_peak_dLight3a6,'tail','right')
[h,p_response_decr_ttest_dLight3a8] = ttest(early_peak_dLight3a8,late_peak_dLight3a8,'tail','right')
[h,p_response_decr_ttest_all] = ttest(early_peak_all,late_peak_all,'tail','right')


% Save data
save(fullfile(filePath, 'early_vs_late_stats.mat'),'p_response_decr_signrank_dLight3a6',...
    'p_response_decr_ttest_dLight3a6','early_peak_dLight3a6_avg',...
    'late_peak_dLight3a6_avg','early_peak_dLight3a6_median','late_peak_dLight3a6_median',...
    'p_response_decr_signrank_dLight3a8','p_response_decr_ttest_dLight3a8','early_peak_dLight3a8_avg',...
    'late_peak_dLight3a8_avg','early_peak_dLight3a8_median','late_peak_dLight3a8_median','early_peak_all_avg',...
    'late_peak_all_avg','p_response_decr_ttest_all');


% Make plot
data_dLight3a6 = [early_peak_dLight3a6 late_peak_dLight3a6];
data_dLight3a8 = [early_peak_dLight3a8 late_peak_dLight3a8];

if horiz_displace == 1
    x_dLight3a6 = [0.4 1.4];
    x_dLight3a8 = [0.6 1.6];
    fig1 = figure;
    hold on
    % Plot data points; left points are early peaks, right points are late
    % peaks
    scatter(x_dLight3a6,data_dLight3a6,'ob','jitter','on','jitterAmount',0.05);
    scatter(x_dLight3a8,data_dLight3a8,'og','jitter','on','jitterAmount',0.05);
    %Plot averages for each sensor
    plot([x_dLight3a6(1)-0.05;x_dLight3a6(1)+0.05],[early_peak_dLight3a6_avg early_peak_dLight3a6_avg],'k-','LineWidth',3)
    plot([x_dLight3a6(2)-0.05;x_dLight3a6(2)+0.05],[late_peak_dLight3a6_avg late_peak_dLight3a6_avg],'k-','LineWidth',3)
    plot([x_dLight3a8(1)-0.05;x_dLight3a8(1)+0.05],[early_peak_dLight3a8_avg early_peak_dLight3a8_avg],'k-','LineWidth',3)
    plot([x_dLight3a8(2)-0.05;x_dLight3a8(2)+0.05],[late_peak_dLight3a8_avg late_peak_dLight3a8_avg],'k-','LineWidth',3)
    
    xlim([0 1.8])
    ylim([min(min(all_data))-0.5 max(max(all_data))+0.5])
    % set(gca,'xtick',[])
    xticklabels({[],[],'early trials',[],[],[],[],'late trials', [],[]})
    % parallelcoords(data_plot,'Labels',labels) %plots lines connecting each pair
    hold off
else
    x = [1 2];
    fig1 = figure;
    hold on
    % Plot data points; left points are early peaks, right points are late
    % peaks
    scatter(x,data_dLight3a6,'b','filled','MarkerEdgeColor','none','jitter','on','jitterAmount',0.05);
    scatter(x,data_dLight3a8,'g','filled','MarkerEdgeColor','none','jitter','on','jitterAmount',0.05);
    % Plot average for all sensors combined
    plot([x(1)-0.1;x(1)+0.1],[early_peak_all_avg early_peak_all_avg],'k-','LineWidth',3)
    plot([x(2)-0.1;x(2)+0.1],[late_peak_all_avg late_peak_all_avg],'k-','LineWidth',3)
    xlim([0 3])
    ylim([min(min(all_data))-0.5 max(max(all_data))+0.5])
    xticklabels({[],[],'early trials',[],'late trials', []})
    hold off
end
figname = 'early vs late comparison scatterplot dLight3';

savefig(fig1,fullfile(filePath,figname))    % save as matlab fig
saveas(fig1,fullfile(filePath,figname),'png')   % save as png for easy import to powerpoint
saveas(fig1,fullfile(filePath,figname),'pdf')   % save as pdf
print(fig1,fullfile(filePath,figname),'-depsc','-tiff')   %saves as eps file


% Plot each dLight version separately
fig2 = figure;
hold on
x = [1 2];
scatter(x,data_dLight3a6,'b','filled','jitter','on','jitterAmount',0.05);
plot([x(1)-0.1;x(1)+0.1],[early_peak_dLight3a6_avg early_peak_dLight3a6_avg],'k-','LineWidth',3)
plot([x(2)-0.1;x(2)+0.1],[late_peak_dLight3a6_avg late_peak_dLight3a6_avg],'k-','LineWidth',3)
xlim([0 3])
dLight3a6_peaks = [early_peak_dLight3a6; late_peak_dLight3a6];
ylim([min(min(all_data))-0.5 max(max(all_data))+0.5])
xticklabels({[],[],'early trials',[],'late trials', []})
hold off
figname = 'early vs late comparison scatterplot dLight3a6';

savefig(fig2,fullfile(filePath,figname))    % save as matlab fig
saveas(fig2,fullfile(filePath,figname),'png')   % save as png for easy import to powerpoint
saveas(fig2,fullfile(filePath,figname),'pdf')   % save as pdf
print(fig2,fullfile(filePath,figname),'-depsc','-tiff')   %saves as eps file


fig3 = figure;
hold on
x = [1 2];
scatter(x,data_dLight3a8,'g','filled','jitter','on','jitterAmount',0.05);
plot([x(1)-0.1;x(1)+0.1],[early_peak_dLight3a8_avg early_peak_dLight3a8_avg],'k-','LineWidth',3)
plot([x(2)-0.1;x(2)+0.1],[late_peak_dLight3a8_avg late_peak_dLight3a8_avg],'k-','LineWidth',3)
xlim([0 3])
dLight3a8_peaks = [early_peak_dLight3a8; late_peak_dLight3a8];
ylim([min(min(all_data))-0.5 max(max(all_data))+0.5])
xticklabels({[],[],'early trials',[],'late trials', []})
hold off
figname = 'early vs late comparison scatterplot dLight3a8';

savefig(fig3,fullfile(filePath,figname))    % save as matlab fig
saveas(fig3,fullfile(filePath,figname),'png')   % save as png for easy import to powerpoint
saveas(fig3,fullfile(filePath,figname),'pdf')   % save as pdf
print(fig3,fullfile(filePath,figname),'-depsc','-tiff')   %saves as eps file






