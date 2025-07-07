close all
clear all

filePath_green = 'X:\2025_05_07\AJ054_SeroSnFR1p4\session03_looming\stack_green';
filePath_red = 'X:\2025_05_07\AJ054_SeroSnFR1p4\session03_looming\stack_red';
identifier = '20250507_03';
file = '20*.tif'; %aligned zstack
nTrials = 80; %80;                                                                                                                   
threshold = 17; %value below which indicates an event
% Look at z profile of red channel to find approx_segment lengths and
% first_trial_indices
approx_segment_lengths = [20 24]; %[20 24]; %[22 22]; %[19 26]; %22; %expected # frames for each segment; determine empirically %protocol change 11/19/22
first_trial_indices = [1 20; 21 44]; %[1  20; 21 43]; %frame index for end of 1st segment of 1st trial; determine empirically

% Set up control for visual stimulus light contamination
filePath_control = 'X:\2025_05_07\AJ054_SeroSnFR1p4\session01_looming_control\stack';
looming_control = 1; %set to 1 if control w/o laser is a looming session
ON_OFF_control = 0; %set to 1 if control w/o laser is a ON/OFF session
nTrials_control = 20; %20; %80 %20
looming_startFrame_control = 19; %19; %frame at which looming period is triggered. Note that actual onset of looming will be 2 frames later.
% NOTE: looming_startFrame_control will usually be 21 for data taken prior
% to 11/19/2022, and 19 for data taken on or after 11/19/2022. Check total
% # frames in session; if 800, looming_startFrame_control=21; if 780,
% looming_startFrame_control=19.

% old_data = 0; %For data taken prior to 11/19/2022 

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

% Create matrix of indices for the first and last frame of each segment of a trial
segmentsPerTrial = length(approx_segment_lengths); % # of segments in each trial
nSegments = nTrials*segmentsPerTrial;
frame_indices = zeros(nSegments,2); %1st col is starting frame of segment, 2nd col is final frame of segment

%%
% Use the expected segment lengths to estimate the frame where a segment
% should end, then look at a few frames around this area to see where the
% red channel signal falls below threshold. Account for the case at the end
% of the session, where the approximate end may exceed the total # of
% frames in the session.
for ii = 1:nTrials
    for jj = 1:segmentsPerTrial
        if ii == 1 %1st trial is special case bcs sensitivity of photodiode is reset
            frame_indices(jj,:) = first_trial_indices(jj,:);
        % The following ugly hack is in cases where segment length is
        % unusual. These are set on a case-by-case basis. Comment out
        % before starting analysis on a new dataset.
        % elseif ii == 34 && jj == 1
        %     frame_indices((ii-1)*segmentsPerTrial+jj,:) = [1455 1469];
        % elseif ii == 60 && jj == 1
        %     frame_indices((ii-1)*segmentsPerTrial+jj,:) = [2588 2610];
        % elseif ii == 47 && jj == 1
        %     frame_indices((ii-1)*segmentsPerTrial+jj,:) = [1934 1952];
        else
            previous_segment_end = frame_indices(ii*segmentsPerTrial+jj-segmentsPerTrial-1,2);
            segment_start = previous_segment_end+1;  
            segment_approx_end = segment_start+[approx_segment_lengths(jj)-3:approx_segment_lengths(jj)+3]; %frame indices close to end of segment 
            if segment_approx_end(end) > size(avg_data_red_chnl,1) && segment_approx_end(end-1) > size(avg_data_red_chnl,1) %deal with case where indices exceed the end of session
                segment_end = segment_approx_end(end-2);
            elseif segment_approx_end(end) > size(avg_data_red_chnl,1) && segment_approx_end(end-1) == size(avg_data_red_chnl,1)
                segment_end = segment_approx_end(end-1);
            else
                segment_end = segment_approx_end(avg_data_red_chnl(segment_approx_end)<threshold)-1; %frame index where next segment starts
            end
            segment_end=segment_end(1); %in case there are >1 indices below threshold, pick the first frame 
            frame_indices((ii-1)*segmentsPerTrial+jj,:) = [segment_start, segment_end];
        end
    end
end

%%
all_segment_lengths = frame_indices(:,2) - frame_indices(:,1) + 1;
min_segment_length = min(all_segment_lengths); %some segments may have extra frames; find the minimum
framesPerTrial = min_segment_length*segmentsPerTrial;

% Create modified version of array "frame_indices" to subtract extra frames
% from segments that are too long
frame_indices_mod = frame_indices; %mod = modified
for jj = 1:nSegments
    this_segment_length = all_segment_lengths(jj);
    if this_segment_length > min_segment_length
        excess = this_segment_length - min_segment_length; % # of extra frames
        frame_ind_thisSegment = frame_indices_mod(jj,:);
        if mod(jj,2) == 1 %odd segment, baseline period
            frame_ind_thisSegment(1) = frame_ind_thisSegment(1) + excess; %add excess to starting index to delete first frame(s)
        elseif mod(jj,2) == 0 %even segment, looming period
            frame_ind_thisSegment(2) = frame_ind_thisSegment(2) - excess; %subtract excess from ending index to delete last frame(s)
        end
        frame_indices_mod(jj,:) = frame_ind_thisSegment;
    end     
end


%% Get fluorescence values from green channel
% Load data from green channel
disp('loading green channel')
data_green_chnl = single(loadImgSequence(filePath_green,file)); %loads zstack as 3d matrix
avg_data_green_chnl = squeeze(mean(mean(data_green_chnl,1),2)); %avg each frame and reduce dimensions

% Apportion the data into segments, using the modified frame indices
% Each row is a segment; each column is a frame
segmented_data = zeros(nSegments,min_segment_length);
for kk = 1:nSegments
    segmented_data(kk,:) = avg_data_green_chnl(frame_indices_mod(kk,1):frame_indices_mod(kk,2))';
end
 


%% Determine stimulus artifact
% The looming protocol begins with a white screen, then a dark circle grows
% until it reaches full-screen size. We much correct for the change in
% luminance, which is detected by the PMTs but should not be included as
% actual dLight signal.

% Get values for raw intensity for control data, averaged across pixels
disp('loading control data')
data_control = single(loadImgSequence(filePath_control,file)); %load control data
data_control_avg = mean(mean(data_control)); %avg all pixel values for each frame
framesTotal_control = size(data_control_avg,3);
framesPerTrial_control = framesTotal_control/nTrials_control;
control_shaped = reshape(data_control_avg,framesPerTrial_control,nTrials_control); %each row is a frame, each col a trial

if looming_control == 1
    % If looming protocol was used as control, divide data into trials
    % where there was a stimulus (odd trials) and trials w/o stim (even
    % trials)
    trialsWithStim_control = [];
    trialsWithNoStim_control = [];
    for jj=1:nTrials_control
        if mod(jj,2)==1 %odd trials
            trialsWithStim_control = [trialsWithStim_control,control_shaped(:,jj)];
        elseif mod(jj,2)==0
            trialsWithNoStim_control = [trialsWithNoStim_control,control_shaped(:,jj)];
        end
    end

    % Average across trials
    trialsWithStim_control_avg = mean(trialsWithStim_control,2); 
    trialsWithNoStim_control_avg = mean(trialsWithNoStim_control,2); 
    

    % Determine the difference btn a completely white screen and completely dark screen; this will be
    % subtracted from the pre-looming period. Each trial with stim starts with
    % white screen and ends with dark screen
    nFramesToAvg = 5; % # of frames to average for white screen and dark screen
    % Take 5 frames from start of trial for white screen. Use frames 2-6
    % instead of 1-5 to avoid any scan initiation artifacts in the first frame.
    white_frames = mean(mean(trialsWithStim_control(2:nFramesToAvg+1,:)));
    % Take 5 frames from end of trial for dark screen. Take frames 6 from
    % the end to one from the end, not the very last frame, in case there
    % is some artifact in the last frame
    dark_frames = mean(mean(trialsWithStim_control(end-nFramesToAvg:end-1,:)));
    % Calculate the difference btn illuminated and dark screen
    artifact_white_screen = abs(white_frames - dark_frames); %diff in luminance btn white and dark screen
   
    % Now determine the artifact during the looming period
    % Before the dark spot reaches full size, fractional decreases in luminance
    % must be considered. nArtifactLoomFrames is the # of frames after
    % looming start during which some fraction of the stimulus artifact
    % occurs. This # must be determined empirically, or made to be long enough to accommodate 
    % almost the entire looming period (e.g. 20 frames).
    nArtifactLoomFrames = floor(framesPerTrial/2)-2; %10;
    artifact_looming = zeros(1,nArtifactLoomFrames); %avg across trials

    for kk = 1:nArtifactLoomFrames-1
        artifact_looming(kk) = abs(trialsWithStim_control_avg(looming_startFrame_control+kk-1,:) - dark_frames);
    end
    
    % Sometimes the first couple values of looming period are artificially
    % high, which appears to be an artifact of scan initiation (since the
    % control protocol triggers acquisition at trial start and looming
    % start). So if these values are larger than the value determined for
    % artifact_white_screen, change them to equal artifact_white_screen.
    artifact_looming(artifact_looming>artifact_white_screen) = artifact_white_screen;
    
    
    

elseif ON_OFF_control == 1
    % If ON/OFF protocol was used as control, take the difference btn light (ON) frames and dark (OFF) frames
    artifact_white_screen = zeros(1,nTrials_control);
    nFramesToAvg = 5; % # of frames to average for white screen and dark screen
    for jj = 1:nTrials_control
%         baseline_frames = mean([mean(control_shaped(1:8,jj)),mean(control_shaped(30:37,jj))]);
%         stimulus_frames = mean(control_shaped(12:25,jj));
        baseline_frames = mean([mean(control_shaped(1:nFramesToAvg,jj)),mean(control_shaped(end-nFramesToAvg:end,jj))]);
        stimulus_frames = mean(control_shaped(round(framesPerTrial_control/2)-nFramesToAvg:round(framesPerTrial_control/2)+nFramesToAvg,jj));
        artifact_white_screen(jj) = abs(stimulus_frames - baseline_frames);
    end
    % Average the difference across trials
    artifact_white_screen = mean(artifact_white_screen);
end


%% Stimulus light artifact must be subtracted from frames where screen was bright (i.e. dark looming circle is absent)
% In array segmented_data, odd rows are the baseline period (white screen); even rows are the looming period; but also
% the even trials are blank, i.e. no looming stimulus during the looming period. Hence every 4th row has the looming 
% stimulus, starting from row 2 (i.e. row 2, 6, 10, 14 etc). For these rows, the stimulus artifact diminishes as the
% looming circle grows to its final size. For all other rows, the stim light artifact should be subtracted from every frame.

% Make copy of unsubtracted data for future reference
% Reshape the data so that # rows = # trials (where segments per trial = 2)
disp('artifact subtraction')
data_raw_unsubtr = zeros(nTrials,framesPerTrial);
for ii = 1:nTrials
    data_raw_unsubtr(ii,1:framesPerTrial/2) = segmented_data(2*ii-1,:);
    data_raw_unsubtr(ii,framesPerTrial/2+1:framesPerTrial) = segmented_data(2*ii,:);
end
% Transpose for consistency with downstream code
data_raw_unsubtr = data_raw_unsubtr'; %now each col is a trial, each row is a frame

% Split into trials with and without looming stim
% Odd # trials have visual stim, even # don't
raw_trialsWithStim_unsubtr = [];
raw_trialsWithNoStim_unsubtr = [];
for jj=1:nTrials
    if mod(jj,2)==1 
        raw_trialsWithStim_unsubtr = [raw_trialsWithStim_unsubtr,data_raw_unsubtr(:,jj)];
    elseif mod(jj,2)==0
        raw_trialsWithNoStim_unsubtr = [raw_trialsWithNoStim_unsubtr,data_raw_unsubtr(:,jj)];  
    end
end

% Average across trials
raw_trialsWithStim_unsubtr_avg = mean(raw_trialsWithStim_unsubtr,2);
raw_trialsWithNoStim_unsubtr_avg = mean(raw_trialsWithNoStim_unsubtr,2);


% Now create another dataset with stimulus light artifact subtracted
% Start with unsubtracted data
raw_trialsWithStim = raw_trialsWithStim_unsubtr;
% Frames prior to looming will be subtracted by artifact_white_screen;
% frames post-looming will be subtracted by looming_artifact. So first we
% need to designate the frame that looming starts. Note that there is a
% delay before looming actually starts bcs it doesn't start until trigger
% offset, but this delay is the same as in the control trials, so the
% timing is still aligned btn experiment and control.
looming_startFrame = min_segment_length+1;
% Subtract artifact_white_screen from pre-looming frames 
raw_trialsWithStim(1:looming_startFrame-1,:) = raw_trialsWithStim(1:looming_startFrame-1,:) - artifact_white_screen;
% If ON_OFF control was used, subtract artifact_white_screen frome first
% frame of looming period. If looming control was used, subtract
% artifact_looming from the appropriate # of frames of looming period.
if ON_OFF_control == 1
    raw_trialsWithStim(framesPerTrial/2+1:framesPerTrial/2+3,:) = raw_trialsWithStim(framesPerTrial/2+1:framesPerTrial/2+3,:) - artifact_white_screen;
%     raw_trialsWithStim(framesPerTrial/2+1,:) = raw_trialsWithStim(framesPerTrial/2+1,:) - artifact_white_screen;
elseif looming_control == 1
    raw_trialsWithStim(looming_startFrame:looming_startFrame+nArtifactLoomFrames-1,:) = raw_trialsWithStim(looming_startFrame:looming_startFrame+nArtifactLoomFrames-1,:) - artifact_looming';
end

% Trials w/o stim get the artifact subtracted from all frames
% Start with unsubtracted data
raw_trialsWithNoStim = raw_trialsWithNoStim_unsubtr;
% Subtract artifact
raw_trialsWithNoStim = raw_trialsWithNoStim - artifact_white_screen;

% Average across trials
raw_trialsWithStim_avg = mean(raw_trialsWithStim,2);
raw_trialsWithNoStim_avg = mean(raw_trialsWithNoStim,2);


%% Calculate dff0
disp('calculating dff0')
method=2; % 1: mode; 2: fixed, e.g. static period; 3: 10% mimimum;
% bck.firstBaselineFrame=framesPerTrial/2-5; %use last 5 trials of baseline (white screen) period
bck.nBaselineFrames=5; %5   %how many frames from each trial to be used as baseline
% bck.framesPerTrial=framesPerTrial; 
bck.avg_intensity_baseline = 0; %set to 1 for avg across pix; not relevant here
% bck.fixPersti=[bck.nBaselineFrames bck.framesPerTrial]; 

% Control trials (stimulus artifact)
% Note that in control protocol, pre-looming and post-looming periods have
% same # frames, and acquisition is triggered, not continuous
bck.firstBaselineFrame=floor(framesPerTrial_control/2)-bck.nBaselineFrames;
bck.framesPerTrial=framesPerTrial_control;
bck.fixPersti=[bck.nBaselineFrames bck.framesPerTrial];
dff0_control = zeros(framesPerTrial_control,nTrials_control/2);
for ii = 1:nTrials_control/2
    dff0_control_thisTrial = calculatedDff0_Kat(trialsWithStim_control(1:framesPerTrial_control,ii),method,bck);
    dff0_control(:,ii) = dff0_control_thisTrial;
end
dff0_control_avg = mean(dff0_control,2);


% For continuous trials, calculate baseline prior to looming onset
bck.framesPerTrial=framesPerTrial;
bck.firstBaselineFrame=framesPerTrial/2-5; %Note: this works bcs trial
% indices are modified such that pre-looming and post-looming segments have
% equal length. However, it may sometimes be better to determine baseline frames with
% respect to looming_startFrame, as shown below. I changed this in August
% 2023 bcs Quartz acquisition rate had slown down.
% bck.firstBaselineFrame=looming_startFrame - bck.nBaselineFrames - 1; %subtract 1 so that baseline frames don't include the frame of looming start, which may be slightly imprecise
% bck.firstBaselineFrame=15;
bck.fixPersti=[bck.nBaselineFrames bck.framesPerTrial];

% Calculate dff0 for data WITHOUT stim artifact subtracted
dff0_trialsWithStim_unsubtr = zeros(framesPerTrial,nTrials/2);
dff0_trialsWithNoStim_unsubtr = zeros(framesPerTrial,nTrials/2);
for ii=1:nTrials/2
    dff0_trialsWithStim_thisTrial = calculatedDff0_Kat(raw_trialsWithStim_unsubtr(:,ii),method,bck);
    dff0_trialsWithStim_unsubtr(:,ii) = dff0_trialsWithStim_thisTrial;
    dff0_trialsWithNoStim_thisTrial = calculatedDff0_Kat(raw_trialsWithNoStim_unsubtr(:,ii),method,bck); 
    dff0_trialsWithNoStim_unsubtr(:,ii) = dff0_trialsWithNoStim_thisTrial;
end

dff0_trialsWithStim_unsubtr_avg = mean(dff0_trialsWithStim_unsubtr,2);
dff0_trialsWithNoStim_unsubtr_avg = mean(dff0_trialsWithNoStim_unsubtr,2);



% Calculate dff0 for data with stim artifact subtracted
dff0_trialsWithStim = zeros(framesPerTrial,nTrials/2);
dff0_trialsWithNoStim = zeros(framesPerTrial,nTrials/2);
for ii=1:nTrials/2
    dff0_trialsWithStim_thisTrial = calculatedDff0_Kat(raw_trialsWithStim(:,ii),method,bck);
    dff0_trialsWithStim(:,ii) = dff0_trialsWithStim_thisTrial;
    dff0_trialsWithNoStim_thisTrial = calculatedDff0_Kat(raw_trialsWithNoStim(:,ii),method,bck);
    dff0_trialsWithNoStim(:,ii) = dff0_trialsWithNoStim_thisTrial;
end

dff0_trialsWithStim_avg = mean(dff0_trialsWithStim,2);
dff0_trialsWithNoStim_avg = mean(dff0_trialsWithNoStim,2);



%%
% Create datasets where outlier trials are omitted
% These are trials where any frame is < or > 3*stdev from the mean
% Outliers will be based on raw data, but will be applied to
% all datasets

% Evaluate each frame of raw_trialsWithStim across trials to find outliers
outliers = [];
for kk = 1:framesPerTrial
    ThisFrame = raw_trialsWithStim(kk,:);
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
raw_trialsWithStim_no_outliers = raw_trialsWithStim;
dff0_trialsWithStim_no_outliers = dff0_trialsWithStim;
raw_trialsWithStim_unsubtr_no_outliers = raw_trialsWithStim_unsubtr;
dff0_trialsWithStim_unsubtr_no_outliers = dff0_trialsWithStim_unsubtr;

% Delete columns (trials) with outliers from all datasets
for n = 1:length(outliers)
    raw_trialsWithStim_no_outliers(:,outliers(n)) = [];
    dff0_trialsWithStim_no_outliers(:,outliers(n)) = [];
    raw_trialsWithStim_unsubtr_no_outliers(:,outliers(n)) = [];
    dff0_trialsWithStim_unsubtr_no_outliers(:,outliers(n)) = [];
end

% Calculate averages
raw_trialsWithStim_no_outliers_avg = mean(raw_trialsWithStim_no_outliers,2);
dff0_trialsWithStim_no_outliers_avg = mean(dff0_trialsWithStim_no_outliers,2);
raw_trialsWithStim_unsubtr_no_outliers_avg = mean(raw_trialsWithStim_unsubtr_no_outliers,2);
dff0_trialsWithStim_unsubtr_no_outliers_avg = mean(dff0_trialsWithStim_unsubtr_no_outliers,2);


% Do the same for trials without stimulus
outliers = [];
for kk = 1:framesPerTrial
    ThisFrame = raw_trialsWithNoStim(kk,:);
    avgThisFrame = mean(ThisFrame);
    stdevThisFrame = std(ThisFrame);
    limit = 3*stdevThisFrame; %values beyond the limit are outliers
    % Find high outlier values
    upper_limit = avgThisFrame + limit;
    high_outlier = find(ThisFrame>upper_limit);
    outliers = [outliers;high_outlier'];
    % Find low outlier values
    lower_limit = avgThisFrame - limit;
    low_outlier = find(ThisFrame<lower_limit);
    outliers = [outliers;low_outlier'];
end
outliers = unique(outliers); %eliminate duplicate values
outliers = sort(outliers,'descend'); %descending order bcs want to delete highest # cols first

% Start with copies of each dataset
raw_trialsWithNoStim_no_outliers = raw_trialsWithNoStim;
dff0_trialsWithNoStim_no_outliers = dff0_trialsWithNoStim;
raw_trialsWithNoStim_unsubtr_no_outliers = raw_trialsWithNoStim_unsubtr;
dff0_trialsWithNoStim_unsubtr_no_outliers = dff0_trialsWithNoStim_unsubtr;

% Delete columns (trials) with outliers from all datasets
for n = 1:length(outliers)
    raw_trialsWithNoStim_no_outliers(:,outliers(n)) = [];
    dff0_trialsWithNoStim_no_outliers(:,outliers(n)) = [];
    raw_trialsWithNoStim_unsubtr_no_outliers(:,outliers(n)) = [];
    dff0_trialsWithNoStim_unsubtr_no_outliers(:,outliers(n)) = [];
end

% Calculate averages
raw_trialsWithNoStim_no_outliers_avg = mean(raw_trialsWithNoStim_no_outliers,2);
dff0_trialsWithNoStim_no_outliers_avg = mean(dff0_trialsWithNoStim_no_outliers,2);
raw_trialsWithNoStim_unsubtr_no_outliers_avg = mean(raw_trialsWithNoStim_unsubtr_no_outliers,2);
dff0_trialsWithNoStim_unsubtr_no_outliers_avg = mean(dff0_trialsWithNoStim_unsubtr_no_outliers,2);






%%
% Make figures

% Plot control with no laser; individual trials with averages
fig1=figure;
hold on

% Looming actually starts 1 sec after the looming period is triggered
% (since the trigger lasts for 1s). Since frame rate is 2.2 Hz, 
% the 1st frame of the looming period is taken 0-0.45s post-trigger, 
% the next frame is taken 0.45-0.9s post-trigger, and looming starts 
% during the acquisition of the 3rd frame.
% So the value of the 3rd frame includes some fraction of looming time. To show
% looming onset for plot, the vertical line will be slightly after the 2nd
% frame, or looming_start_control+1.
% Distance between x axis points (frames) is 0.45s, and looming starts 0.1s
% after 2nd frame: 0.1s/0.45s = 0.2. 
% So looming onset is looming_startFrame_control+1+0.2
looming_start_control = looming_startFrame_control+1.2;

if looming_control == 1
    for ii=1:nTrials_control/2
        plot(trialsWithNoStim_control(:,ii),'Color',[0.5 0.5 0.5])
        p1 = plot(trialsWithStim_control(:,ii),'-r');
        p1.Color = [p1.Color 0.5]; %make lines semi-transparent
    end
    plot(trialsWithNoStim_control_avg,'k','LineWidth',3)
    plot(trialsWithStim_control_avg,'r','LineWidth',3)
    yLim = get(gca,'yLim');  
    plot([looming_start_control looming_start_control],yLim,'b--')
    
elseif ON_OFF_control == 1
    for ii=1:nTrials_control
        p1 = plot(control_shaped(:,ii),'-r');
        p1.Color = [p1.Color 0.5]; %make lines semi-transparent
    end
    raw_fluor_avg = mean(control_shaped,2);
    plot(raw_fluor_avg,'r','LineWidth',3)
    yLim = get(gca,'yLim');  
    plot([10 10],yLim,'b--')
    plot([29 29],yLim,'b--')
end
    
hold off
figname = 'raw fluor no laser control';
save_figures(fig1,filePath_green,result,figname)


%% Data with stim artifact subtracted
% Plot individual trials with averages, raw fluorescence


fig2=figure;
hold on

for ii=1:nTrials/2
    plot(raw_trialsWithNoStim(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(raw_trialsWithStim(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(raw_trialsWithNoStim_avg,'k','LineWidth',3)
plot(raw_trialsWithStim_avg,'r','LineWidth',3)

yLim = get(gca,'yLim'); 

% See note in plot 1 about how to designate looming onset.
looming_start = looming_startFrame+1.2;
plot([looming_start looming_start],yLim,'b--') %vert line when looming starts

hold off
figname = 'raw fluor with and without stimulus';
save_figures(fig2,filePath_green,result,figname)


% Plot individual trials with averages, dff0
fig3=figure;
hold on

for ii=1:nTrials/2    
    plot(dff0_trialsWithNoStim(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(dff0_trialsWithStim(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(dff0_trialsWithNoStim_avg,'k','LineWidth',3)
plot(dff0_trialsWithStim_avg,'r','LineWidth',3)

yLim = get(gca,'yLim');  
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'dff0 with and without stimulus';
save_figures(fig3,filePath_green,result,figname)



% Plot averages with standard error
% Calculate standard error for each frame, across trials
stError_trialsWithStim = zeros(1,framesPerTrial);
stError_trialsWithNoStim = zeros(1,framesPerTrial);

for k = 1:framesPerTrial 
    stError_trialsWithStim(k) = nanstd(dff0_trialsWithStim(k,:))/sqrt(nTrials/2);
    stError_trialsWithNoStim(k) = nanstd(dff0_trialsWithNoStim(k,:))/sqrt(nTrials/2);
end

fig4 = figure;
hold on
props1 = {'-k','lineWidth',2};
props2 = {'-r','lineWidth',3};
shadedErrorBar([],dff0_trialsWithNoStim_avg,stError_trialsWithNoStim,props1); %empty brackets for x axis
shadedErrorBar([],dff0_trialsWithStim_avg,stError_trialsWithStim,props2); %empty brackets for x axis

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'avg sterror with and wo stim';
save_figures(fig4,filePath_green,result,figname)


% Plot the difference between trials with and without visual stim
stimNoStim_difference = dff0_trialsWithStim_avg - dff0_trialsWithNoStim_avg;
fig5 = figure;
hold on
plot(stimNoStim_difference,'b','LineWidth',3)

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'difference with and wo stim';
save_figures(fig5,filePath_green,result,figname)
    
    

%% Data WITHOUT stim artifact subtracted
% Plot individual trials with averages, raw fluorescence
fig6=figure;
hold on

for ii=1:nTrials/2
    plot(raw_trialsWithNoStim_unsubtr(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(raw_trialsWithStim_unsubtr(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(raw_trialsWithNoStim_unsubtr_avg,'k','LineWidth',3)
plot(raw_trialsWithStim_unsubtr_avg,'r','LineWidth',3)

yLim = get(gca,'yLim'); 
plot([looming_start looming_start],yLim,'b--') %vert line when looming starts

hold off
figname = 'raw fluor with and without stimulus unsubtr';
save_figures(fig6,filePath_green,result,figname)


% Plot individual trials with averages, dff0
fig7=figure;
hold on

for ii=1:nTrials/2    
    plot(dff0_trialsWithNoStim_unsubtr(:,ii),'Color',[0.5 0.5 0.5])
    p1 = plot(dff0_trialsWithStim_unsubtr(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(dff0_trialsWithNoStim_unsubtr_avg,'k','LineWidth',3)
plot(dff0_trialsWithStim_unsubtr_avg,'r','LineWidth',3)

yLim = get(gca,'yLim');  
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'dff0 with and without stimulus unsubtr';
save_figures(fig7,filePath_green,result,figname)



% Plot averages with standard error
% Calculate standard error for each frame, across trials
stError_trialsWithStim_unsubtr = zeros(1,framesPerTrial);
stError_trialsWithNoStim_unsubtr = zeros(1,framesPerTrial);

for k = 1:framesPerTrial 
    stError_trialsWithStim_unsubtr(k) = nanstd(dff0_trialsWithStim_unsubtr(k,:))/sqrt(nTrials/2);
    stError_trialsWithNoStim_unsubtr(k) = nanstd(dff0_trialsWithNoStim_unsubtr(k,:))/sqrt(nTrials/2);
end

fig8 = figure;
hold on
props1 = {'-k','lineWidth',2};
props2 = {'-r','lineWidth',3};
shadedErrorBar([],dff0_trialsWithNoStim_unsubtr_avg,stError_trialsWithNoStim_unsubtr,props1); %empty brackets for x axis
shadedErrorBar([],dff0_trialsWithStim_unsubtr_avg,stError_trialsWithStim_unsubtr,props2); %empty brackets for x axis

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'avg sterror with and wo stim unsubtr';
save_figures(fig8,filePath_green,result,figname)

% Plot the difference between trials with and without visual stim
stimNoStim_difference_unsubtr = dff0_trialsWithStim_unsubtr_avg - dff0_trialsWithNoStim_unsubtr_avg;
fig9 = figure;
hold on
plot(stimNoStim_difference_unsubtr,'b','LineWidth',3)

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'difference with and wo stim unsubtr';
save_figures(fig9,filePath_green,result,figname)


%% Figures for datasets where outliers are omitted
% Use datasets where stim artifact is subtracted

% Plot individual trials with averages, raw fluorescence
fig10=figure;
hold on

for ii=1:size(raw_trialsWithNoStim_no_outliers,2)
    plot(raw_trialsWithNoStim_no_outliers(:,ii),'Color',[0.5 0.5 0.5])
end
for ii=1:size(raw_trialsWithStim_no_outliers,2)
    p1 = plot(raw_trialsWithStim_no_outliers(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(raw_trialsWithNoStim_no_outliers_avg,'k','LineWidth',3)
plot(raw_trialsWithStim_no_outliers_avg,'r','LineWidth',3)

yLim = get(gca,'yLim'); 
plot([looming_start looming_start],yLim,'b--') %vert line when looming starts

hold off
figname = 'raw fluor with and without stimulus no outliers';
save_figures(fig10,filePath_green,result,figname)


% Plot individual trials with averages, dff0
fig11=figure;
hold on

for ii=1:size(dff0_trialsWithNoStim_no_outliers,2)    
    plot(dff0_trialsWithNoStim_no_outliers(:,ii),'Color',[0.5 0.5 0.5])
end
for ii=1:size(dff0_trialsWithStim_no_outliers,2)    
    p1 = plot(dff0_trialsWithStim_no_outliers(:,ii),'-r');
    p1.Color = [p1.Color 0.5]; %make lines semi-transparent
end
plot(dff0_trialsWithNoStim_no_outliers_avg,'k','LineWidth',3)
plot(dff0_trialsWithStim_no_outliers_avg,'r','LineWidth',3)

yLim = get(gca,'yLim');  
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'dff0 with and without stimulus no outliers';
save_figures(fig11,filePath_green,result,figname)


% Plot averages with standard error
% Calculate standard error for each frame, across trials
stError_trialsWithStim_no_outliers = zeros(1,framesPerTrial);
stError_trialsWithNoStim_no_outliers = zeros(1,framesPerTrial);

for k = 1:framesPerTrial 
    stError_trialsWithStim_no_outliers(k) = nanstd(dff0_trialsWithStim_no_outliers(k,:))/sqrt(size(dff0_trialsWithStim_no_outliers,2)/2);
    stError_trialsWithNoStim_no_outliers(k) = nanstd(dff0_trialsWithNoStim_no_outliers(k,:))/sqrt(size(dff0_trialsWithStim_no_outliers,2)/2);
end

fig12 = figure;
hold on
props1 = {'-k','lineWidth',2};
props2 = {'-r','lineWidth',3};
shadedErrorBar([],dff0_trialsWithNoStim_no_outliers_avg,stError_trialsWithNoStim_no_outliers,props1); %empty brackets for x axis
shadedErrorBar([],dff0_trialsWithStim_no_outliers_avg,stError_trialsWithStim_no_outliers,props2); %empty brackets for x axis

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'avg sterror with and wo stim no outliers';
save_figures(fig12,filePath_green,result,figname)


% Plot the difference between trials with and without visual stim
stimNoStim_difference_no_outliers = dff0_trialsWithStim_no_outliers_avg - dff0_trialsWithNoStim_no_outliers_avg;
fig13 = figure;
hold on
plot(stimNoStim_difference_no_outliers,'b','LineWidth',3)

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'difference with and wo stim no outliers';
save_figures(fig13,filePath_green,result,figname)


% Plot raw fluorescence for control and unsubtracted dLight data on same plot
% To account for the difference in baseline fluorescence btn control and
% dLight, avg the first 10 frames (blank screen) of each trace and subtract
% this from all values of each trace. This sets the "baseline" (pre-looming
% values) to zero.
control_baseline_avg = mean(trialsWithStim_control_avg(1:10));
dLight_baseline_avg = mean(raw_trialsWithStim_unsubtr_no_outliers_avg(1:10));
control_baseline_subtr = trialsWithStim_control_avg - control_baseline_avg;
dLight_baseline_subtr = raw_trialsWithStim_unsubtr_no_outliers_avg - dLight_baseline_avg;

% Find standard error
stError_raw_control = zeros(1,size(trialsWithStim_control_avg,1));
for k = 1:size(trialsWithStim_control_avg,1) 
    stError_raw_control(k) = nanstd(trialsWithStim_control(k,:))/sqrt(size(trialsWithStim_control,2)/2);
end

stError_raw_trialsWithStim_unsubtr = zeros(1,size(raw_trialsWithStim_unsubtr_no_outliers_avg,1));
for k = 1:size(raw_trialsWithStim_unsubtr_no_outliers_avg,1) 
    stError_raw_trialsWithStim_unsubtr(k) = nanstd(raw_trialsWithStim_unsubtr(k,:))/sqrt(size(raw_trialsWithStim_unsubtr,2)/2);
end

% Align the traces to 1st frame of looming period. Since # frames for control can be different than # frames for dLight,
% chop off the extra frames from either the beginning of the trial
% (pre-looming period) or end of the trial (looming period)
% Set up vecytors to be plotted, which will be shortened depending on
% differences in # frames:
dLight_plot = dLight_baseline_subtr;
control_plot = control_baseline_subtr;
stError_dLight_plot = stError_raw_trialsWithStim_unsubtr;
stError_control_plot = stError_raw_control;
% Look at pre-looming period
dLight_preloom = dLight_baseline_subtr(1:looming_startFrame-1);
control_preloom = control_baseline_subtr(1:looming_startFrame_control-1);
extra_frames_preloom = abs(length(dLight_preloom) - length(control_preloom));
if length(dLight_preloom) > length(control_preloom)
    dLight_plot = dLight_baseline_subtr(1+extra_frames_preloom:end);
    stError_dLight_plot = stError_raw_trialsWithStim_unsubtr(1+extra_frames_preloom:end);
else
    control_plot = control_baseline_subtr(1+extra_frames_preloom:end);
    stError_control_plot = stError_raw_control(1+extra_frames_preloom:end);
end
% Look at post-looming period
dLight_postloom = dLight_baseline_subtr(looming_startFrame:end);
control_postloom = control_baseline_subtr(looming_startFrame_control:end);
extra_frames_postloom = abs(length(dLight_postloom) - length(control_postloom));
if length(dLight_postloom) > length(control_postloom)
    dLight_plot = dLight_plot(1:end-extra_frames_postloom);
    stError_dLight_plot = stError_dLight_plot(1:end-extra_frames_postloom);
else
    control_plot = control_plot(1:end-extra_frames_postloom);
    stError_control_plot = stError_control_plot(1:end-extra_frames_postloom);
end
% We also need to adjust the point where looming onset is plotted.
% Onset will be the next frame after the length of the shortest preloom
% period, plus 1.2 (see explanation for fig1)
if length(dLight_preloom) > length(control_preloom)
    looming_start_adj = length(control_preloom)+1+1.2;
else
    looming_start_adj = length(dLight_preloom)+1+1.2;
end


% Make the figure
fig14 = figure;
hold on

props1 = {'-k','lineWidth',2};
props2 = {'-r','lineWidth',3};
shadedErrorBar([],dLight_plot,stError_dLight_plot,props2); %empty brackets for x axis
shadedErrorBar([],control_plot,stError_control_plot,props1); %empty brackets for x axis
% Re-plot the main line for dLight data for visibility
plot(dLight_plot,'r','LineWidth',3)

yLim = get(gca,'yLim');
plot([looming_start_adj looming_start_adj],yLim,'b--')

hold off
figname = 'avg raw unsubtracted with stim vs control';
save_figures(fig14,filePath_green,result,figname)


% Plot difference in raw fluorescence btn dLight and control
dLight_control_diff = dLight_plot-control_plot;
fig15 = figure;
hold on
plot(dLight_control_diff,'b','LineWidth',3)

yLim = get(gca,'yLim');
plot([looming_start looming_start],yLim,'b--')

hold off
figname = 'difference raw fluor dLight vs control';
save_figures(fig15,filePath_green,result,figname)



%% Save data
save(fullfile(filePath_green, ['analysis_' identifier '.mat']),'all_segment_lengths','approx_segment_lengths',...
    'artifact_white_screen','bck','dff0_trialsWithNoStim','dff0_trialsWithNoStim_avg',...
    'dff0_trialsWithNoStim_no_outliers','dff0_trialsWithNoStim_no_outliers_avg','dff0_trialsWithNoStim_unsubtr',...
    'dff0_trialsWithNoStim_unsubtr_avg','dff0_trialsWithNoStim_unsubtr_no_outliers','dff0_trialsWithNoStim_unsubtr_no_outliers_avg',...
    'dff0_trialsWithStim','dff0_trialsWithStim_avg','dff0_trialsWithStim_no_outliers','dff0_trialsWithStim_no_outliers_avg',...
    'dff0_trialsWithStim_unsubtr','dff0_trialsWithStim_unsubtr_avg','dff0_trialsWithStim_unsubtr_no_outliers',...
    'dff0_trialsWithStim_unsubtr_no_outliers_avg','frame_indices','frame_indices_mod','framesPerTrial',...
    'framesPerTrial_control','looming_control','looming_start','looming_startFrame','nTrials','nTrials_control','ON_OFF_control',...
    'raw_trialsWithNoStim','raw_trialsWithNoStim_avg','raw_trialsWithNoStim_no_outliers',...
    'raw_trialsWithNoStim_no_outliers_avg','raw_trialsWithNoStim_unsubtr','raw_trialsWithNoStim_unsubtr_avg',...
    'raw_trialsWithNoStim_unsubtr_no_outliers','raw_trialsWithNoStim_unsubtr_no_outliers_avg','raw_trialsWithStim',...
    'raw_trialsWithStim_avg','raw_trialsWithStim_no_outliers','raw_trialsWithStim_no_outliers_avg','raw_trialsWithStim_unsubtr',...
    'raw_trialsWithStim_unsubtr_avg','raw_trialsWithStim_unsubtr_no_outliers','raw_trialsWithStim_unsubtr_no_outliers_avg',...
    'stError_trialsWithStim','stError_trialsWithNoStim','stError_trialsWithStim_unsubtr','stError_trialsWithNoStim_unsubtr',...
    'stError_trialsWithStim_no_outliers','stError_trialsWithNoStim_no_outliers','threshold')





%%
function save_figures(fig,filePath,result,figname)
savefig(fig,fullfile(filePath,result,figname))    % save as matlab fig
saveas(fig,fullfile(filePath,result,figname),'png')   % save as png for easy import to powerpoint
% print(fig,fullfile(filePath,figname),'-depsc','-tiff')   %saves as eps file
end




