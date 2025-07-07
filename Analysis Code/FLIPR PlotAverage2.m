close all
clear all

relative_data = 1; %if 1 -> data is calculated relative

% mice = { 'TVK347L' 'TVK348L' 'TVK349L' 'TVK350L' 'TVK347R' 'TVK348R' 'TVK349R' 'TVK350R'};
% Power_used = [1 0.5 1 0.6 1 0.8 0.8 0.8];
% AFscalingfactor = [1 0.25 1 0.364 1 0.62 0.62 0.62];

mice = { 'TVK348L' 'TVK349L' 'TVK350L' 'TVK348R' 'TVK349R' 'TVK350R'};
Power_used = [ 0.5 1 0.6 0.8 0.8 0.8];
AFscalingfactor = [ 0.25 1 0.364  0.62 0.62 0.62];

traceNames = {};
for j = 1:length(mice)
    mouseName = mice{j}
    load([mice{j} '_data'])
    mouseLtData = tempSave;

    daqAlignmentIdx = mouseLtData.daqAlignmentIdx;
    Lifetime = mouseLtData.lifetime;
    DC = mouseLtData.VDC;
    mDC(j) = mean(DC);

    downSampledTo = mouseLtData.settings.downSampledTo;
 
    timewindow = 10;
    
    for z = 1:length(daqAlignmentIdx)
        if daqAlignmentIdx(z) < 1
            daqAlignmentIdx = daqAlignmentIdx(1:z-1);
        end
    end

    %Align DC
    for jj = 1:length(daqAlignmentIdx)
        trialDataDC(jj,:) = DC( daqAlignmentIdx(jj) - timewindow * downSampledTo : daqAlignmentIdx(jj) + timewindow*downSampledTo);
    end
    
    %Calculate df/f
    tBaseline = 5;
    trialDataDC = trialDataDC *-1;
    for q = 1:size(trialDataDC,1)
        F = mean(trialDataDC(q,1:tBaseline*downSampledTo));
        trialDataDCdf(q,:) = (trialDataDC(q,:) - F)  ./ F;
    end

    %plot DC
    figure()
    plot(trialDataDCdf')
    hold on
    
    title(['Fluorescence intensity individual trials ' mouseLtData.mouse])
    ylabel('Fluorescence lifetime(ns)')
    xlabel('Time')
    xticks([0:1*downSampledTo:20*downSampledTo])
    xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'})

    sessionAverageDC = mean(trialDataDCdf,1);
    figure(100)
    hold on
    plot(sessionAverageDC, 'DisplayName', mouseLtData.mouse)

    title('fluorescence intensity ')
    ylabel('Fluorescence intensity (df/f)')
    xlabel('Time')
    xticks([0:1*downSampledTo:20*downSampledTo])
    % xticklabels({'-5' '-4' '-3' '-2' '-1' '0' '1' '2' '3' '4' '5'})
    xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'}) 
    

    %Plot lifetimes
    for jj = 1:length(daqAlignmentIdx)
        trialData(jj,:) = Lifetime( daqAlignmentIdx(jj) - timewindow * downSampledTo : daqAlignmentIdx(jj) + timewindow*downSampledTo);
    end

    meanBaseline = mean(trialData(:,1:timewindow/2*downSampledTo), 'all');
    
    for jj = 1:size(trialData,1)
        baselinePerTrial = mean(trialData(jj,1:timewindow/2*downSampledTo), 2);
        trialNorm(jj,:) = (trialData(jj,:) - baselinePerTrial)*1000;
        maxPeak(j,jj) = max(trialNorm(jj,180:220));
    end
    
    figure()
    plot(trialData'*1000 - meanBaseline*1000)
    hold on
    
    title(['Fluorescence lifetime individual trials ' mouseLtData.mouse])
    ylabel('Fluorescence lifetime(ns)')
    xlabel('Time')
    xticks([0:1*downSampledTo:20*downSampledTo])
    xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'})

    sessionAverage = mean(trialData,1);
    figure(1000)
    hold on
    plot(sessionAverage, 'DisplayName', mouseLtData.mouse)

    title('fluorescence lifetime ')
    ylabel('Fluorescence lifetime(ns)')
    xlabel('Time')
    xticks([0:1*downSampledTo:20*downSampledTo])
    % xticklabels({'-5' '-4' '-3' '-2' '-1' '0' '1' '2' '3' '4' '5'})
    xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'})

    mouseLtData.trialData = trialData;
    mouseLtData.timewindow = timewindow;
    mouseLtData.sessionAverage = sessionAverage;


    if relative_data == 1
        deltaLt = sessionAverage - meanBaseline;
    else
        deltaLt = sessionAverage;
    end

    figure(1001)
    hold on
    plot(deltaLt, 'DisplayName', mouseLtData.mouse)

    title('fluorescence lifetime ')
    ylabel('Fluorescence lifetime(ns)')
    xlabel('Time')
    xticks([0:1*downSampledTo:20*downSampledTo])
    % xticklabels({'-5' '-4' '-3' '-2' '-1' '0' '1' '2' '3' '4' '5'})
    xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'})

    mouseLtData.trialData = trialData;
    mouseLtData.timewindow = timewindow;
    mouseLtData.sessionAverage = sessionAverage;

    save([mice{j} '_data.mat'], "tempSave")
    
    allMice.cohortMiceSessionAverage(j,:) = sessionAverage;
    allMice.cohortMiceDeltaLt(j,:) = deltaLt;
    allMice.cohortMiceVDC(j,:) = mean(mouseLtData.VDC);
    allMice.sessionAverageDC(j,:) = sessionAverageDC;
    
end


%% create average trace

avrALlMice = mean(allMice.cohortMiceDeltaLt*1000);
stdAllMice = std(allMice.cohortMiceDeltaLt*1000,1);

figure(10001)
y = avrALlMice; % your mean vector;
x = 1:numel(y);
std_dev = stdAllMice;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.7 0.7 0.7]);
hold on;
plot(x, y, 'k', 'LineWidth', 2);
hold on

xticks([0:1*downSampledTo:20*downSampledTo])
xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'})
title('Average delta lifetime across mice')
ylabel('Fluorescence lifetime(ps)')
xlabel('Time(s)')


avrALlMiceDC = mean(allMice.sessionAverageDC);
stdAllMiceDC = std(allMice.sessionAverageDC,1);

figure(10002)
y = avrALlMiceDC; % your mean vector;
x = 1:numel(y);
std_dev = stdAllMiceDC;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.7 0.7 0.7]);
hold on;
plot(x, y, 'k', 'LineWidth', 2);
hold on

xticks([0:1*downSampledTo:20*downSampledTo])
xticklabels({'-10' '' '-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10'})
title('Average delta intensity across mice')
ylabel('Fluorescence intensity(df/f)')
xlabel('Time(s)')





