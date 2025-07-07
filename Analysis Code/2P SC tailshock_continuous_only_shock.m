close all
clear all

%%Structure of stimulus: 
% 40 trials; odd trials have shock, even trials are blank


filePath_green = 'X:\2025_05_07\AJ054_SeroSnFR1p4\session05_tailshock\stack_green';
filePath_red = 'X:\2025_05_07\AJ054_SeroSnFR1p4\session05_tailshock\stack_red';
identifier = '20250507_05';
file = '20*.tif'; %aligned zstack
nTrials = 40; %total trials in entire session
threshold = 15; %value below which indicates an event

approx_segment_lengths = [22,66]; %[22,66]; %expected # frames for each segment; determine empirically
first_trial_indices = [1  22; 23 88]; %[1 22; 23 88]; %frame index for end of 1st segment of 1st trial; determine empirically
% Create the result directory if it doesn't exist
result = 'plots full frame';
if ~exist(fullfile(filePath_green,result))
    mkdir(fullfile(filePath_green,result));
end


%% Figure out frame indices of each segment of the trials
% Load data from red channel
disp('loading red channel')
data_red_chnl = single(loadImgSequence(filePath_red,file)); %loads zstack as 3d matrix
avg_data_red_chnl = squeeze(mean(mean(data_red_chnl,1),2)); %avg each frame and reduce dimensions
nFrames = length(avg_data_red_chnl);


% Create matrix of indices for the first and last frame of each segment of a trial, 
segmentsPerTrial = length(approx_segment_lengths); % # of segments in each trial
nSegments = nTrials*segmentsPerTrial;
frame_indices = zeros(nSegments,2); %1st col is starting frame of segment, 2nd col is final frame of segment

% Find frame indices for trials 
for ii = 1:nTrials
    for jj = 1:segmentsPerTrial
        if ii == 1 %1st trial is special case bcs sensitivity of photodiode is reset
            frame_indices(jj,:) = first_trial_indices(jj,:);
        % elseif ii == 37 && jj==2
        %     frame_indices((ii-1)*segmentsPerTrial+jj,:) = [3085 3148];
        % elseif ii == 38 && jj==1
        %     frame_indices((ii-1)*segmentsPerTrial+jj,:) = [3148 3168];
        else
            previous_segment_end = frame_indices(ii*segmentsPerTrial+jj-segmentsPerTrial-1,2);
            segment_start = previous_segment_end+1; 
            segment_approx_end = [segment_start+approx_segment_lengths(jj)-3,segment_start+approx_segment_lengths(jj)-2, segment_start+approx_segment_lengths(jj)-1, segment_start+approx_segment_lengths(jj), segment_start+approx_segment_lengths(jj)+1, segment_start+approx_segment_lengths(jj)+2,  segment_start+approx_segment_lengths(jj)+3]; %frame indices close to end of segment  
            segment_end = segment_approx_end(find(avg_data_red_chnl(segment_approx_end)<threshold))-1; %frame index where next segment starts
            segment_end=segment_end(1); %in case there are >1 indices below threshold, pick the first frame 
            frame_indices((ii-1)*segmentsPerTrial+jj,:) = [segment_start, segment_end];
        end
    end
end

% Since the # of frames may differ from trial to trial, need to correct the
% frame indices
all_segment_lengths = frame_indices(:,2) - frame_indices(:,1) + 1;
all_segment_lengths = reshape(all_segment_lengths,segmentsPerTrial,[]);
min_segment_lengths = min(all_segment_lengths,[],2); %min of each row to get shortest length of each segment

% Create modified version of frame indices array to subtract extra frames
% from segments that are too long
frame_indices_mod = zeros(nSegments,2); %mod = modified

for mm = 1:nTrials
    for kk = 1:segmentsPerTrial
        this_segment_length = all_segment_lengths(kk,mm);
        if this_segment_length > min_segment_lengths(kk)
            excess = this_segment_length -  min_segment_lengths(kk);
            frame_ind_thisSegment = frame_indices((mm-1)*segmentsPerTrial+kk,:);
            if mod(kk,2) == 1 %odd segment, presetim period
                frame_ind_thisSegment(1) = frame_ind_thisSegment(1) + excess; %add 1 to starting index to delete first frame
            elseif mod(kk,2) == 0 %even segment, stim period or poststim period
                frame_ind_thisSegment(2) = frame_ind_thisSegment(2) - excess; %subtract 1 from ending index to delete last frame
            end
            frame_indices_mod((mm-1)*segmentsPerTrial+kk,:) = frame_ind_thisSegment;
        else
            frame_indices_mod((mm-1)*segmentsPerTrial+kk,:) = frame_indices((mm-1)*segmentsPerTrial+kk,:);
        end     
    end
end

framesPerTrial = floor(frame_indices_mod(end)/nSegments);

%% Get fluorescence values from green channel
% Load data from green channel
disp('loading green channel')
data_green_chnl = single(loadImgSequence(filePath_green,file)); %loads zstack as 3d matrix
avg_data_green_chnl = squeeze(mean(mean(data_green_chnl,1),2)); %avg each frame and reduce dimensions

% Apportion the data into segments, using the modified frame indices
% Each row is a trial; each column is a frame
disp('analyzing trials')
segmented_data = zeros(nTrials,sum(min_segment_lengths));
for mm = 1:nTrials
    indicesThisTrial = frame_indices_mod((mm-1)*segmentsPerTrial+1:(mm-1)*segmentsPerTrial+segmentsPerTrial,:);
    dataThisTrial = [];
    for kk = 1:segmentsPerTrial
        dataThisSegment = avg_data_green_chnl(indicesThisTrial(kk,1):indicesThisTrial(kk,2));
        dataThisTrial = [dataThisTrial;dataThisSegment];
    end
    segmented_data(mm,:) = dataThisTrial;
end



% Split into trials with and w/o stimulus
% Array segmented_data: each row is a trial; cols are frames
raw_trialsWithShock = [];
raw_trialsWithNoShock = [];
for mm = 1:nTrials
    if mod(mm,2)==1 %odd trials have stim
        raw_trialsWithShock = [raw_trialsWithShock;segmented_data(mm,:)];
    elseif mod(mm,2)==0 %even trials don't have stim
        raw_trialsWithNoShock = [raw_trialsWithNoShock;segmented_data(mm,:)];  
    end
end

% Transpose for consistency with downstream code (dff0 calculation)
raw_trialsWithShock = raw_trialsWithShock'; %now each col is a trial, each row is a frame
raw_trialsWithNoShock = raw_trialsWithNoShock';
% Average across trials
raw_trialsWithShock_avg = mean(raw_trialsWithShock,2);
raw_trialsWithNoShock_avg = mean(raw_trialsWithNoShock,2);


%% Calculate dff0
method=2; % 1: mode; 2: fixed, e.g. static period; 3: 10% mimimum;
bck.nBaselineFrames=5; %5   %how many frames from each trial to be used as baseline
bck.firstBaselineFrame=min_segment_lengths(1) - bck.nBaselineFrames - 1; %use last few frames of prestim segment for baseline; subtract 1 to avoid the very last frame, 
%since timing of frame acquisition relative to stimulus onset can be slightly variable in continuous recording.
bck.framesPerTrial=sum(min_segment_lengths); 
bck.avg_intensity_baseline = 0; %set to 1 for avg across pix; not relevant here
bck.fixPersti=[bck.nBaselineFrames bck.framesPerTrial]; 

% Trials with vis+aud stim and shock
dff0_trialsWithShock = zeros(size(raw_trialsWithShock,1),size(raw_trialsWithShock,2));
dff0_trialsWithNoShock = zeros(size(raw_trialsWithShock,1),size(raw_trialsWithNoShock,2));
for ii=1:nTrials/2
    dff0_trialsWithShock_thisTrial = calculatedDff0_Kat(raw_trialsWithShock(:,ii),method,bck);
    dff0_trialsWithShock(:,ii) = dff0_trialsWithShock_thisTrial;
    dff0_trialsWithNoShock_thisTrial = calculatedDff0_Kat(raw_trialsWithNoShock(:,ii),method,bck);
    dff0_trialsWithNoShock(:,ii) = dff0_trialsWithNoShock_thisTrial;
end
dff0_trialsWithShock_avg = mean(dff0_trialsWithShock,2);
dff0_trialsWithNoShock_avg = mean(dff0_trialsWithNoShock,2);


%%
% Create datasets where outlier trials are omitted
% These are trials where any frame is < or > 3*stdev from the mean
% Outliers will be based on raw data, but will be applied to
% all datasets

% Evaluate each frame of raw_trialsWithShock across trials to find outliers
outliers = [];
for kk = 1:framesPerTrial
    ThisFrame = raw_trialsWithShock(kk,:);
    avgThisFrame = mean(ThisFrame);
    stdevThisFrame = std(ThisFrame);
    limit = 3*stdevThisFrame; %values beyond the limit are outliers
    % Find high outlier values
    upper_limit = avgThisFrame + limit;
    high_outlier = find(ThisFrame>upper_limit);
    outliers = [outliers;high_outlier']; %need ' in case there's >1 trials
    % Find low outlier values
    lower_limit = avgThisFrame - limit;
    low_outlier = find(ThisFrame<lower_limit);
    outliers = [outliers;low_outlier'];
end
outliers = unique(outliers); %eliminate duplicate values
outliers = sort(outliers,'descend'); %descending order bcs want to delete highest # cols first

% Start with copies of each dataset
raw_trialsWithShock_no_outliers = raw_trialsWithShock;
dff0_trialsWithShock_no_outliers = dff0_trialsWithShock;
% Delete columns (trials) with outliers from all datasets
for n = 1:length(outliers)
    raw_trialsWithShock_no_outliers(:,outliers(n)) = [];
    dff0_trialsWithShock_no_outliers(:,outliers(n)) = [];
end
% Calculate averages
raw_trialsWithShock_no_outliers_avg = mean(raw_trialsWithShock_no_outliers,2);
dff0_trialsWithShock_no_outliers_avg = mean(dff0_trialsWithShock_no_outliers,2);


% Do the same for trials without shock
outliers = [];
for kk = 1:framesPerTrial
    ThisFrame = raw_trialsWithNoShock(kk,:);
    avgThisFrame = mean(ThisFrame);
    stdevThisFrame = std(ThisFrame);
    limit = 3*stdevThisFrame; %values beyond the limit are outliers
    % Find high outlier values
    upper_limit = avgThisFrame + limit;
    high_outlier = find(ThisFrame>upper_limit);
    outliers = [outliers;high_outlier']; %need ' in case there's >1 trials
    % Find low outlier values
    lower_limit = avgThisFrame - limit;
    low_outlier = find(ThisFrame<lower_limit);
    outliers = [outliers;low_outlier'];
end
outliers = unique(outliers); %eliminate duplicate values
outliers = sort(outliers,'descend'); %descending order bcs want to delete highest # cols first

% Start with copies of each dataset
raw_trialsWithNoShock_no_outliers = raw_trialsWithNoShock;
dff0_trialsWithNoShock_no_outliers = dff0_trialsWithNoShock;
% Delete columns (trials) with outliers from all datasets
for n = 1:length(outliers)
    raw_trialsWithNoShock_no_outliers(:,outliers(n)) = [];
    dff0_trialsWithNoShock_no_outliers(:,outliers(n)) = [];
end
% Calculate averages
raw_trialsWithNoShock_no_outliers_avg = mean(raw_trialsWithNoShock_no_outliers,2);
dff0_trialsWithNoShock_no_outliers_avg = mean(dff0_trialsWithNoShock_no_outliers,2);



%% Save data
save(fullfile(filePath_green, ['analysis_' identifier '.mat']),'raw_trialsWithShock','raw_trialsWithNoShock',...
'bck','frame_indices','frame_indices_mod','dff0_trialsWithShock','dff0_trialsWithNoShock',...
'framesPerTrial','all_segment_lengths','min_segment_lengths','segmented_data',...
'raw_trialsWithShock_no_outliers','raw_trialsWithNoShock_no_outliers',...
'dff0_trialsWithShock_no_outliers','dff0_trialsWithNoShock_no_outliers','raw_trialsWithShock_avg',...
'raw_trialsWithNoShock_avg','dff0_trialsWithShock_avg','dff0_trialsWithNoShock_avg',...
'raw_trialsWithShock_no_outliers_avg','raw_trialsWithNoShock_no_outliers_avg',...
'dff0_trialsWithShock_no_outliers_avg','dff0_trialsWithNoShock_no_outliers_avg');


%%
% Make figures

% Plot individual trials with averages, raw fluorescence
fig1=figure;
hold on

for ii=1:nTrials/2
    plot(raw_trialsWithNoShock(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(raw_trialsWithShock(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(raw_trialsWithNoShock_avg,'k','LineWidth',3)
plot(raw_trialsWithShock_avg,'r','LineWidth',3)

yLim = get(gca,'yLim'); 
shock_start = min_segment_lengths(1)+1;
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'raw fluorescence';
savefig(fig1,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig1,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


% Plot individual trials with averages, dff0
fig2=figure;
hold on

for ii=1:nTrials/2    
    plot(dff0_trialsWithNoShock(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(dff0_trialsWithShock(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(dff0_trialsWithNoShock_avg,'k','LineWidth',3)
plot(dff0_trialsWithShock_avg,'r','LineWidth',3)

yLim = get(gca,'yLim');  
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'dff0';
savefig(fig2,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig2,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


% Plot averages with standard error
% Calculate standard error for each frame, across trials
stError_trialsWithShock = zeros(1,sum(min_segment_lengths));
stError_trialsWithNoShock = zeros(1,sum(min_segment_lengths));

for k = 1:sum(min_segment_lengths) 
    stError_trialsWithShock(k) = nanstd(dff0_trialsWithShock(k,:))/sqrt(nTrials/2);
    stError_trialsWithNoShock(k) = nanstd(dff0_trialsWithNoShock(k,:))/sqrt(nTrials/2);
end

fig3 = figure;
hold on
props1 = {'-k','lineWidth',2};
props2 = {'-r','lineWidth',3};
shadedErrorBar([],dff0_trialsWithNoShock_avg,stError_trialsWithNoShock,props1); %empty brackets for x axis
shadedErrorBar([],dff0_trialsWithShock_avg,stError_trialsWithShock,props2); %empty brackets for x axis

yLim = get(gca,'yLim');
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'avg sterror';
savefig(fig3,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig3,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


% Plot the difference between trials with and without visual stim
shockNOshock_difference = dff0_trialsWithShock_avg - dff0_trialsWithNoShock_avg;
fig4 = figure;
hold on
plot(shockNOshock_difference,'k','LineWidth',3)

yLim = get(gca,'yLim');
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'difference with and wo shock';
savefig(fig4,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig4,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


% Plot data with outlier trials removed, raw fluorescence
fig5 = figure;
hold on
nTrials_no_outliers = size(raw_trialsWithShock_no_outliers,2);

for ii=1:nTrials_no_outliers
    plot(raw_trialsWithNoShock_no_outliers(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(raw_trialsWithShock_no_outliers(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(raw_trialsWithNoShock_no_outliers_avg,'k','LineWidth',3)
plot(raw_trialsWithShock_no_outliers_avg,'r','LineWidth',3)

yLim = get(gca,'yLim'); 
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'raw fluorescence no outliers';
savefig(fig5,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig5,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint  


% Plot data with outlier trials removed, dff0
fig6 = figure;
hold on
dff0_trialsWithNoShock_no_outliers_avg = mean(dff0_trialsWithNoShock_no_outliers,2);
dff0_trialsWithShock_no_outliers_avg = mean(dff0_trialsWithShock_no_outliers,2);


for ii=1:nTrials_no_outliers
    plot(dff0_trialsWithNoShock_no_outliers(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(dff0_trialsWithShock_no_outliers(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(dff0_trialsWithNoShock_no_outliers_avg,'k','LineWidth',3)
plot(dff0_trialsWithShock_no_outliers_avg,'r','LineWidth',3)

yLim = get(gca,'yLim'); 
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'dff0 no outliers';
savefig(fig6,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig6,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


% Plot averages with standard error, outliers removed
stError_trialsWithShock_no_outliers = zeros(1,sum(min_segment_lengths));
stError_trialsWithNoShock_no_outliers = zeros(1,sum(min_segment_lengths));

for k = 1:sum(min_segment_lengths) 
    stError_trialsWithShock_no_outliers(k) = nanstd(dff0_trialsWithShock_no_outliers(k,:))/sqrt(nTrials_no_outliers/2);
    stError_trialsWithNoShock_no_outliers(k) = nanstd(dff0_trialsWithNoShock_no_outliers(k,:))/sqrt(nTrials_no_outliers/2);
end

fig7 = figure;
hold on
props1 = {'-k','lineWidth',2};
props2 = {'-r','lineWidth',3};
shadedErrorBar([],dff0_trialsWithNoShock_no_outliers_avg,stError_trialsWithNoShock_no_outliers,props1); %empty brackets for x axis
shadedErrorBar([],dff0_trialsWithShock_no_outliers_avg,stError_trialsWithShock_no_outliers,props2); %empty brackets for x axis

yLim = get(gca,'yLim');
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'avg sterror no outliers';
savefig(fig7,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig7,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


% Plot the difference between trials with and without visaud stim, outliers
% removedload
shockNOshock_difference_no_outliers = dff0_trialsWithShock_no_outliers_avg - dff0_trialsWithNoShock_no_outliers_avg;
fig8 = figure;
hold on
plot(shockNOshock_difference_no_outliers,'k','LineWidth',3)

yLim = get(gca,'yLim');
plot([shock_start shock_start],yLim,'b--')

hold off
figname = 'difference with and wo shock no outliers';
savefig(fig8,fullfile(filePath_green,result,figname))    % save as matlab fig
saveas(fig8,fullfile(filePath_green,result,figname),'png')   % save as png for easy import to powerpoint


