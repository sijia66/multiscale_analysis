% major changes to the code
% instead of using all trials, just use trials
% missing scripts: calcSeqGuessCount, calcSeqGuess, calcSeqReachDuration
%   
% 

%m181022_basicSeqProperties
clear all

global MONKEYDIR
MONKEYDIR = 'C:/Users/Orsborn Lab/OneDrive - UW/projects/Brain EEG/data';


%data properties
days_randOnly = {'180329', '180401', '180408','180415'}; %{'180721', '180722', '180723', '180724', '180726', '180727', '180728', '180729', ...
days_seq = {};
%days_seq = {'180330', '180331', '180401', '180403', '180404', '180405'};
DAYS = [days_randOnly days_seq];
seqDayIDs = find(ismember(DAYS, days_seq));

%load up all trials across days
TRIALS = [];
for d = 1:length(DAYS)
    
    base = [MONKEYDIR '/' DAYS{d} '/mat/'];
    load([base 'Trials.mat'], 'Trials') %just look at Trials
    
    %fn = fieldnames(Trials);
    TRIALS = cat(2, TRIALS, Trials);
end
clear Trials

nTr = length(TRIALS);


%get basic properties of trials
seqRand = [TRIALS.SequenceRandom];
seqLen  = [TRIALS.SequenceLength];
seqLearning = [TRIALS.SequenceLearning];
sTr     = [TRIALS.Success];
%%

%get basic RT, reach time etc. metrics across trials
guessCount = calcSeqGuessCount(TRIALS); % # gueses
[guessPos, guessTargEst, guessTime, guessTimeOff, guessTrID, guessSegID] = calcSeqGuess(TRIALS); %properties of guesses

[acqRel2Cue_times, acqRel2Cue_targids] = calcSeqReachAcqRel2Cue(TRIALS); %timing of acquisition rel. to cue (to see if correctly getting target)

[reachDur_times, reachDur_targids] = calcSeqReachDuration(TRIALS); %duration of reach dur (hand transport)

[acqPos, acq_targids] = calcSeqReachPos(TRIALS); %touch end-points

[reachRT_times, reachRT_targids] = calcSeqReachRT(TRIALS);

[trDur] = calcSeqTrialDur(TRIALS);

[TargPos, TargColor, CenterPos, CenterColor] = calcSeqTargArray(TRIALS);

[cueDelay, targids] = calcSeqCueDelay(TRIALS);


%fix all the target positions b/c labview trial alignment is tricky and not
%working properly yet.
%use calcSeqClosestTarget to fix. Work within-day in case target array
%moves across days.
trDay= {TRIALS.Day};
dayID = nan(nTr,1);

TargPosEst = nan(size(acqPos));
for d=1:length(DAYS)
    inds = ismember(trDay, DAYS{d});
    dayID(inds) = d;
    
    
    TargPosEst(inds,:,:) = calcSeqClosestTarget(acqPos(inds,:,:), acq_targids(inds,:,:), TargPos(inds,:,:), CenterPos(inds,:));
end

%now fix ones where sequence isn't random by forcing to array position
for d=1:length(DAYS)
    inds = dayID'==d & seqLearning==1 & seqRand==0;
    dayID(inds) = d;
   
    TargPosEst(inds,:,:) = calcSeqClosestTarget_predictablePattern(acq_targids(inds,:,:), TargPos(inds,:,:), CenterPos(inds,:));
end

%compute target errors for each position
acqErr = nan(size(acqPos));
for tr=1:nTr
    
    for t=1:size(acqPos,2)
        %         if mod(t,2)==1
        %             acqErr(tr,t,:) = sq(acqPos(tr,t,:)) - CenterPos(tr,:)';
        %         else
        %             acqErr(tr,t,:) = acqPos(tr,t,:) - TargPos(tr,t,:);
        %         end
        acqErr(tr,t,:) = acqPos(tr,t,:) - TargPosEst(tr,t,:);
    end
end
acqErrMag = sqrt( sum( acqErr.^2,3));


%also get indication of whether day is random-only or sequence training day
randDayID = nan(nTr,1);
seqDayID  = nan(nTr,1);

for d=1:length(days_randOnly)
    inds = ismember(trDay, days_randOnly{d});
    randDayID(inds) = d;
end

for d=1:length(days_seq)
    inds = ismember(trDay, days_seq{d});
    seqDayID(inds) = d;
end



X = nanmean(cueDelay(:,1:2:end),2); %average cue delay across targets
cuedTrial = X<300;

randInds = seqRand==1 & seqLearning==1 & seqLen == 2 & sTr==1;
seqInds_len3 = seqRand==0 & seqLearning==1 & seqLen == 2 & ismember(dayID', seqDayIDs) & sTr==1;


%% plot some things to take a look at progressions
% 
% 
% [axH, figH] = plotDataByDay(trDur(randInds), dayID(randInds), 40);
% title(axH, 'Random-only trial duration')
% 
% 
% [axH, figH] = plotDataByDay(trDur(seqInds_len3), dayID(seqInds_len3), 40);
% title(axH, 'Seq len = 3, trial duration')
% 
% 
% [mDur_rand, sdDur_rand, day_list_rand] = calcMetricByDay(trDur(randInds), dayID(randInds), 'median', 'all');
% [mDur_seq, sdDur_seq, day_list_seq] = calcMetricByDay(trDur(seqInds_len3), dayID(seqInds_len3), 'mean', 50, 'first');
% 
% %
% [axH, figH] = plotDataByDay(acqRel2Cue_times(seqInds_len3,2), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(acqRel2Cue_times(seqInds_len3,4), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(acqRel2Cue_times(seqInds_len3,6), dayID(seqInds_len3), 40);
% 
% 
% 
% 
% %
% [axH, figH] = plotDataByDay(reachDur_times(seqInds_len3,1), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(reachDur_times(seqInds_len3,3), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(reachDur_times(seqInds_len3,5), dayID(seqInds_len3), 40);
% 
% 
% 
% [axH, figH] = plotDataByDay(guessCount(seqInds_len3,2), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(guessCount(seqInds_len3,4), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(guessCount(seqInds_len3,6), dayID(seqInds_len3), 40);
% 
% 
% 
% [axH, figH] = plotDataByDay(acqErrMag(seqInds_len3,2), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(acqErrMag(seqInds_len3,4), dayID(seqInds_len3), 40);
% 
% 
% 
% [axH, figH] = plotDataByDay(reachRT_times(seqInds_len3,1), dayID(seqInds_len3), 40);
% [axH, figH] = plotDataByDay(reachRT_times(seqInds_len3,3), dayID(seqInds_len3), 40);
% 
% 
% acqBeforeCue = acqRel2Cue_times < -50;
% 
% [axH, figH] = plotDataByDay(acqBeforeCue(seqInds_len3,2), dayID(seqInds_len3), 40, 1);
% [axH, figH] = plotDataByDay(acqBeforeCue(seqInds_len3,4), dayID(seqInds_len3), 40, 1);
% 
% 
% %success rate
% seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen == 3 & ismember(dayID', seqDayIDs);
% [axH, figH] = plotDataByDay(sTr(seqInds_len3_all)', dayID(seqInds_len3_all), 40, 1, [], [], [0 1]);
% title(axH, 'Seq len = 3, success rate')
% 
% randInds_all = seqRand==1 & seqLearning==1 & seqLen == 3;
% [axH, figH] = plotDataByDay(sTr(randInds_all)', dayID(randInds_all), 40, 1, [], [], [0 1]);
% title(axH, 'Random-only success rate')


%% start making more careful plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
winSize = 20;
fontSize = 14;

colors = cbrewer('seq', 'GnBu', max(dayID)+2);
colors = colors(3:end,:);

%% success rates
figH = figure;
ax1 = subplot(2,5, [1:3]);
ax2 = subplot(2,5, 4);
ax3 = subplot(2,5, 5);
ax4 = subplot(2,5, [6:8]);
ax5 = subplot(2,5, 9);
ax6 = subplot(2,5, 10);


% %plot across-day trajectories for rand only and seq in long axes (1 and 4)
randInds_all = seqRand==1 & seqLearning==1 & seqLen == 3;
% plotDataByDay(sTr(randInds_all)', dayID(randInds_all), winSize, 1, ax1, colors, [0 1]);
% title(ax1, 'Random-only success rate')
% 
% seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen == 2 & ismember(dayID', seqDayIDs);
% plotDataByDay(sTr(seqInds_len3_all)', dayID(seqInds_len3_all), winSize, 1, ax4, colors, [0 1]);
% title(ax4, 'Seq len = 2, success rate')


%plot session-average across days next
[M, SD, x] = calcMetricByDay(sTr(randInds_all)', dayID(randInds_all), 'bino', 'all');
errorbar(ax2, x, M, M-SD, M+SD, '-o','lineWidth', 3);
set(ax2, 'ylim', [0 1], 'xlim', [x(1)-.5, x(end)+.5])

[M, SD, x] = calcMetricByDay(sTr(seqInds_len3_all)', dayID(seqInds_len3_all), 'bino', 'all');
errorbar(ax5, x, M, M-SD, M+SD, '-o', 'lineWidth', 3);
set(ax5, 'ylim', [0 1], 'xlim', [x(1)-.5, x(end)+.5])

%now plot late on each day vs early on next day (to look for savings)
[M, SD, x] = calcMetricByDay_earlyLate(sTr(randInds_all)', dayID(randInds_all), 'bino', 50);
Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
Y(2,:)   = M(1,2:end);   %early for days 2:N
SD_(1,:) = SD(2,1:end-1);
SD_(2,:) = SD(1,2:end);
figure(figH)
set(figH,'currentAxes', ax3)
hold on
for d=1:size(Y,2)
    errorbar([1 2], Y(:,d), Y(:,d)-SD_(:,d), Y(:,d)+SD_(:,d), '-o', 'color', colors(d,:))
end
set(ax3, 'ylim', [0 1], 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})


[M, SD, x] = calcMetricByDay_earlyLate(sTr(seqInds_len3_all)', dayID(seqInds_len3_all), 'bino', 50);
Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
Y(2,:)   = M(1,2:end);   %early for days 2:N
SD_(1,:) = SD(2,1:end-1);
SD_(2,:) = SD(1,2:end);
figure(figH)
set(figH,'currentAxes', ax6)
hold on
for d=1:size(Y,2)
    errorbar([1 2], Y(:,d), Y(:,d)-SD_(:,d), Y(:,d)+SD_(:,d), '-o', 'color', colors(d,:))
end
set(ax6, 'ylim', [0 1], 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})



%% trial duration
figH = figure;
ax1 = subplot(2,5, [1:3]);
ax2 = subplot(2,5, 4);
ax3 = subplot(2,5, 5);
ax4 = subplot(2,5, [6:8]);
ax5 = subplot(2,5, 9);
ax6 = subplot(2,5, 10);

DATA = trDur;
yaxlim = [2000 5000];

winSize = 20;
%inds of data to include
randInds_all = seqRand==1 & seqLearning==1 & seqLen == 2 & sTr==1;
seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen == 2 & ismember(dayID', seqDayIDs) & sTr==1 & ~cuedTrial';


%plot across-day trajectories for rand only and seq in long axes (1 and 4)
plotDataByDay(DATA(randInds_all), dayID(randInds_all), winSize, 0, ax1, colors, yaxlim);
title(ax1, 'Random-only trial duration')

plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 0, ax4, colors, yaxlim);
title(ax4, 'Seq len = 3, trial duration')


%plot session-average across days next
[M, SD, x] = calcMetricByDay(DATA(randInds_all)', dayID(randInds_all), 'mean', 'all');
errorbar(ax2, x, M, SD, SD, '-o','lineWidth', 3);
set(ax2, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])

[M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 'all');
errorbar(ax5, x, M, SD, SD, '-o', 'lineWidth', 3);
set(ax5, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])

%now plot late on each day vs early on next day (to look for savings)
[M, SD, x] = calcMetricByDay_earlyLate(DATA(randInds_all)', dayID(randInds_all), 'mean', 20);
Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
Y(2,:)   = M(1,2:end);   %early for days 2:N
SD_(1,:) = SD(2,1:end-1);
SD_(2,:) = SD(1,2:end);
figure(figH)
set(figH,'currentAxes', ax3)
hold on
for d=1:size(Y,2)
    errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
end
set(ax3, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})


[M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 20);
Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
Y(2,:)   = M(1,2:end);   %early for days 2:N
SD_(1,:) = SD(2,1:end-1);
SD_(2,:) = SD(1,2:end);
figure(figH)
set(figH,'currentAxes', ax6)
hold on
for d=1:size(Y,2)
    errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
end
set(ax6, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})


%% acq rel2cue
targInds = [2 4];
nTarg = length(targInds);

figH = figure;
for t=1:nTarg
    
    eval(['ax' num2str((t-1)*3+1) ' = subplot(nTarg,5, [1:3]+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+2) ' = subplot(nTarg,5, 4+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+3) ' = subplot(nTarg,5, 5+(t-1)*5);'])
    
    
    
    DATA = acqRel2Cue_times(:,targInds(t));
    DATA(DATA<-50) = -50;
    yaxlim = [-75 400];
    
    winSize = 20;
    %inds of data to include
    seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen >= t & ismember(dayID', seqDayIDs) & sTr==1 & ~cuedTrial';
    
    eval(['ax = ax' num2str((t-1)*3+1) ';']);
    plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 0, ax, colors, yaxlim);
    title(ax, ['Seq len = 3, acqisition relative to cue T=' num2str(t)])
    
    
    %plot session-average across days next
    eval(['ax = ax' num2str((t-1)*3+2) ';']);
    [M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 'all');
    errorbar(ax, x, M, SD, SD, 'k-o', 'lineWidth', 3);
    set(ax, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])
    
    %now plot late on each day vs early on next day (to look for savings)
    eval(['ax = ax' num2str((t-1)*3+3) ';']);
    [M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 20);
    Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
    Y(2,:)   = M(1,2:end);   %early for days 2:N
    SD_(1,:) = SD(2,1:end-1);
    SD_(2,:) = SD(1,2:end);
    figure(figH)
    set(figH,'currentAxes', ax)
    hold on
    for d=1:size(Y,2)
        errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
    end
    set(ax, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})
    
end %targets


%% reach duration
targInds = [1 3];
nTarg = length(targInds);

figH = figure;
for t=1:nTarg
    
    eval(['ax' num2str((t-1)*3+1) ' = subplot(nTarg,5, [1:3]+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+2) ' = subplot(nTarg,5, 4+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+3) ' = subplot(nTarg,5, 5+(t-1)*5);'])
    
    
    
    DATA = reachDur_times(:,targInds(t));
    %DATA(DATA<-50) = -50;
    yaxlim = [100 700];
    
    winSize = 20;
    %inds of data to include
    seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen >= t & ismember(dayID', seqDayIDs) & sTr==1 & ~cuedTrial';
    noGuess=guessCount(:,targInds(t)+1)==0;
    seqInds_len3_all = seqInds_len3_all & noGuess';
    
    eval(['ax = ax' num2str((t-1)*3+1) ';']);
    plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 0, ax, colors, yaxlim);
    title(ax, ['Seq len = 3, reach time T=' num2str(t)])
    
    
    %plot session-average across days next
    eval(['ax = ax' num2str((t-1)*3+2) ';']);
    [M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 'all');
    errorbar(ax, x, M, SD, SD, 'k-o', 'lineWidth', 3);
    set(ax, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])
    
    %now plot late on each day vs early on next day (to look for savings)
    eval(['ax = ax' num2str((t-1)*3+3) ';']);
    [M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 20);
    Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
    Y(2,:)   = M(1,2:end);   %early for days 2:N
    SD_(1,:) = SD(2,1:end-1);
    SD_(2,:) = SD(1,2:end);
    figure(figH)
    set(figH,'currentAxes', ax)
    hold on
    for d=1:size(Y,2)
        errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
    end
    set(ax, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})
    
end %targets


%% reach duration, random only
targInds = [1 3];
nTarg = length(targInds);

figH = figure;
for t=1:nTarg
    
    eval(['ax' num2str((t-1)*3+1) ' = subplot(nTarg,5, [1:3]+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+2) ' = subplot(nTarg,5, 4+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+3) ' = subplot(nTarg,5, 5+(t-1)*5);'])
    
    
    
    DATA = reachDur_times(:,targInds(t));
    %DATA(DATA<-50) = -50;
    yaxlim = [100 500];
    
    winSize = 20;
    %inds of data to include
    seqInds_len3_all = seqRand==1 & seqLearning==1 & seqLen == 3 & ismember(dayID', seqDayIDs) & sTr==1;
%     noGuess=guessCount(:,targInds(t)+1)==0;
%     seqInds_len3_all = seqInds_len3_all & noGuess';
%     
    eval(['ax = ax' num2str((t-1)*3+1) ';']);
    plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 0, ax, colors, yaxlim);
    title(ax, ['random, reach time T=' num2str(t)])
    
    
    %plot session-average across days next
    eval(['ax = ax' num2str((t-1)*3+2) ';']);
    [M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 'all');
    errorbar(ax, x, M, SD, SD, 'k-o', 'lineWidth', 3);
    set(ax, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])
    
    %now plot late on each day vs early on next day (to look for savings)
    eval(['ax = ax' num2str((t-1)*3+3) ';']);
    [M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 20);
    Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
    Y(2,:)   = M(1,2:end);   %early for days 2:N
    SD_(1,:) = SD(2,1:end-1);
    SD_(2,:) = SD(1,2:end);
    figure(figH)
    set(figH,'currentAxes', ax)
    hold on
    for d=1:size(Y,2)
        errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
    end
    set(ax, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})
    
end %targets

%% reach error 
targInds = [2 4];
nTarg = length(targInds);


acqErrMag_tmpFix = acqErrMag;
acqErrMag_tmpFix(acqErrMag>10) = nan;
figH = figure;
for t=1:nTarg
    
    eval(['ax' num2str((t-1)*3+1) ' = subplot(nTarg,5, [1:3]+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+2) ' = subplot(nTarg,5, 4+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+3) ' = subplot(nTarg,5, 5+(t-1)*5);'])
    
    
    
    DATA = acqErrMag_tmpFix(:,targInds(t));
    %DATA(DATA<-50) = -50;
    yaxlim = [.5 2.5];
    
    winSize = 20;
    %inds of data to include
    seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen >= t & ismember(dayID', seqDayIDs) & sTr==1 & ~cuedTrial';
    noGuess=guessCount(:,targInds(t))==0;
    seqInds_len3_all = seqInds_len3_all & noGuess';
    
    eval(['ax = ax' num2str((t-1)*3+1) ';']);
    plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 0, ax, colors, yaxlim);
    title(ax, ['Seq len = 3, reach time T=' num2str(t)])
    
    
    %plot session-average across days next
    eval(['ax = ax' num2str((t-1)*3+2) ';']);
    [M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 'all');
    errorbar(ax, x, M, SD, SD, 'k-o', 'lineWidth', 3);
    set(ax, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])
    
    %now plot late on each day vs early on next day (to look for savings)
    eval(['ax = ax' num2str((t-1)*3+3) ';']);
    [M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 20);
    Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
    Y(2,:)   = M(1,2:end);   %early for days 2:N
    SD_(1,:) = SD(2,1:end-1);
    SD_(2,:) = SD(1,2:end);
    figure(figH)
    set(figH,'currentAxes', ax)
    hold on
    for d=1:size(Y,2)
        errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
    end
    set(ax, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})
    
end %targets

%% reach error, random only 
targInds = [2 4];
nTarg = length(targInds);


acqErrMag_tmpFix = acqErrMag;
acqErrMag_tmpFix(acqErrMag>10) = nan;
figH = figure;
for t=1:nTarg
    
    eval(['ax' num2str((t-1)*3+1) ' = subplot(nTarg,5, [1:3]+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+2) ' = subplot(nTarg,5, 4+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+3) ' = subplot(nTarg,5, 5+(t-1)*5);'])
    
    
    
    DATA = acqErrMag_tmpFix(:,targInds(t));
    %DATA(DATA<-50) = -50;
    yaxlim = [.5 20];
    
    winSize = 20;
    %inds of data to include
    seqInds_len3_all = seqRand==1 & seqLearning==1 & seqLen >=t & ismember(dayID', seqDayIDs) & sTr==1;
    noGuess=guessCount(:,targInds(t))==0;
    seqInds_len3_all = seqInds_len3_all & noGuess';
    
    eval(['ax = ax' num2str((t-1)*3+1) ';']);
    plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 0, ax, colors, yaxlim);
    title(ax, ['Seq len = 3, reach time T=' num2str(t)])
    
    
    %plot session-average across days next
    eval(['ax = ax' num2str((t-1)*3+2) ';']);
    [M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 'all');
    errorbar(ax, x, M, SD, SD, 'k-o', 'lineWidth', 3);
    set(ax, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])
    
    %now plot late on each day vs early on next day (to look for savings)
    eval(['ax = ax' num2str((t-1)*3+3) ';']);
    [M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'mean', 20);
    Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
    Y(2,:)   = M(1,2:end);   %early for days 2:N
    SD_(1,:) = SD(2,1:end-1);
    SD_(2,:) = SD(1,2:end);
    figure(figH)
    set(figH,'currentAxes', ax)
    hold on
    for d=1:size(Y,2)
        errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
    end
    set(ax, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})
    
end %targets


%% acq before cue
targInds = [2 4];
nTarg = length(targInds);

acqBeforeCue = acqRel2Cue_times <= -50;

figH = figure;
for t=1:nTarg
    
    eval(['ax' num2str((t-1)*3+1) ' = subplot(nTarg,5, [1:3]+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+2) ' = subplot(nTarg,5, 4+(t-1)*5);'])
    eval(['ax' num2str((t-1)*3+3) ' = subplot(nTarg,5, 5+(t-1)*5);'])
    
    
    
    DATA = acqBeforeCue(:,targInds(t));
    %DATA(DATA<-50) = -50;
    yaxlim = [0 1];
    
    winSize = 20;
    %inds of data to include
    seqInds_len3_all = seqRand==0 & seqLearning==1 & seqLen >= t & ismember(dayID', seqDayIDs) & sTr==1 & ~cuedTrial';
    %noGuess=guessCount(:,targInds(t)+1)==0;
    %seqInds_len3_all = seqInds_len3_all & noGuess';
    
    eval(['ax = ax' num2str((t-1)*3+1) ';']);
    plotDataByDay(DATA(seqInds_len3_all), dayID(seqInds_len3_all), winSize, 1, ax, colors, yaxlim);
    title(ax, ['Seq len = 3, reach time T=' num2str(t)])
    
    
    %plot session-average across days next
    eval(['ax = ax' num2str((t-1)*3+2) ';']);
    [M, SD, x] = calcMetricByDay(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'bino', 'all');
    errorbar(ax, x, M, SD, SD, 'k-o', 'lineWidth', 3);
    set(ax, 'ylim', yaxlim, 'xlim', [x(1)-.5, x(end)+.5])
    
    %now plot late on each day vs early on next day (to look for savings)
    eval(['ax = ax' num2str((t-1)*3+3) ';']);
    [M, SD, x] = calcMetricByDay_earlyLate(DATA(seqInds_len3_all)', dayID(seqInds_len3_all), 'bino', 20);
    Y(1,:)   = M(2,1:end-1); %late for days 1:N-1
    Y(2,:)   = M(1,2:end);   %early for days 2:N
    SD_(1,:) = SD(2,1:end-1);
    SD_(2,:) = SD(1,2:end);
    figure(figH)
    set(figH,'currentAxes', ax)
    hold on
    for d=1:size(Y,2)
        errorbar([1 2], Y(:,d), SD_(:,d), SD_(:,d), '-o', 'color', colors(d,:))
    end
    set(ax, 'ylim', yaxlim, 'xlim', [.5 2.5], 'xtick', [1 2], 'xticklabels', {'L, dN', 'E, dN+1'})
    
end %targets