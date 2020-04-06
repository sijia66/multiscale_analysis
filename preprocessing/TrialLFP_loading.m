
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global MONKEYDIR
MONKEYDIR = 'C:/Users/Orsborn Lab/OneDrive - UW/projects/Brain EEG/data';

%load \180328\mat\trials
%load an experiment file from the same day for loading the electrode
%locations etc. 

% 1 for random sequence trial, this option will create a new
% Trial structure and overwrite the old one
% 0 for all trials
trialInfo.sequence_random = 1; 
trialInfo.filterTime = 0; % filter out the trials that have Reach Time beyond the 75th percentile
trialInfo.filterAcq = 1;
trialInfo.filterAcq_time = 5     ; % second filter out any trial that time elapse between reachStop and TargetAcq
trialInfo.loadFirst100trials = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if trialInfo.sequence_random
    v = [Trials.SequenceRandom];
    Trials = Trials(v ==1);
    
    disp('Only selected  trials with random sequence targets!')
end
    
if trialInfo.filterTime
    startToStop = [Trials.ReachStop] - [Trials.ReachStart];
    
    startToStop_lower  = prctile(startToStop,25);
    startToStop_upper = prctile(startToStop,75);
    
    v = (startToStop <= startToStop_upper) ;
    
    Trials = Trials(v);
    disp('Only selected  trials with lower 75% of the target reach times!')
end

if trialInfo.filterAcq
    stopToAcq = [Trials.ReachStop] - [Trials.TargAq];
    v = abs(stopToAcq) <= trialInfo.filterAcq_time;
    Trials = Trials(v);
    disp('Only selected trials with 5 s between ReachStop to TargAcq!')
end

userinput = input('Press enter to load the trials in the next code section!');

if trialInfo.loadFirst100trials
Trials = Trials(1:100);
disp('Will load the first 100 trials')

end

%%
%inspect the electrodes

%%
trialInfo.tBefore = 0;
trialInfo.tAfter  = 2e2;
trialInfo.trig = '';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;

%load the experiment from the recording
P_chamber_SC32 = ...
    drmap_coordTransform_microdrive2chamber(experiment, 'LM1_SC32_1');
P_chamber_SC32 = P_chamber_SC32';

P_chamber_ECOG_3 = ...
    drmap_coordTransform_microdrive2chamber(experiment, 'LM1_ECOG_3');
P_chamber_ECOG_3 = P_chamber_ECOG_3';

P_chamber = P_chamber_ECOG_3;
%load two drives
driveNameAnalyze = 'LM1_ECOG_3';
drive_idx = 1; %1 for ECOG and 2 for SC32

dataParam.eList = [experiment.hardware.microdrive(drive_idx).electrodes.channelid];
data = trialLfp(Trials,driveNameAnalyze,dataParam.eList,[],trialInfo.trig,[trialInfo.tBefore,trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);
data_twoDrives{1} = data;

driveNameAnalyze = 'LM1_SC32_1';
drive_idx = 2; %1 for ECOG and 2 for SC32
dataParam.eList = [experiment.hardware.microdrive(drive_idx).electrodes.channelid];
data = trialLfp(Trials,driveNameAnalyze,dataParam.eList,[],trialInfo.trig,[trialInfo.tBefore,trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);

data_twoDrives{2} = data;

%% sort out the targets based on positions

[Pos, targids] = calcSeqReachTarg(Trials);
%so far, we can only look athe first activity

Pos_Seq_1 = squeeze(Pos(:,2,:));

%load target format

[filename,path] = uigetfile("","Select target format to load");
selectedfile = fullfile(path,filename);

if path == 0
    
    Pos_Seq_1_unique = unique(Pos_Seq_1,'rows');
    disp('Use target format as they appear in the data')
else
    load(selectedfile)
    disp('loaded selected target format')
end


if size(Pos_Seq_1_unique,1) > 8
    Pos_Seq_1_unique = Pos_Seq_1_unique(1:2:end,:);
    disp('selected every other target')
end 

subplot(1,2,1)
scatter(Pos_Seq_1_unique(:,1),Pos_Seq_1_unique(:,2));
xlabel('mm')
ylabel('mm')
axis equal
hold on
for ti = 1:size(Pos_Seq_1_unique,1)
    text(Pos_Seq_1_unique(ti,1)+0.5,Pos_Seq_1_unique(ti,2),num2str(ti),'FontSize',20);
end
hold off

subplot(1,2,2)
% assignTargetNumbers
trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
histogram(trial_TargetAssigned,7) 
xlabel('Target Number')
ylabel('Counts')


%% 2 sort trials by targets, okay 
%Sorting the tasks! 
dataTarget = mulSortTrials(data_ref_mean,trial_TargetAssigned, Trials);
disp('Sorted grid averaged data in the same structure')

dataTarget_SC32 = mulSortTrials(data_twoDrives{2},trial_TargetAssigned, Trials);

%dataTarget_matrix_min = mulSortTrials_matrix_min(dataTarget);

%% concatenate all data from the trials
driveNameAnalyze = 'LM1_SC32_1';
data_all_trials = concat_all_trials(Trials,experiment,driveNameAnalyze,MONKEYDIR);

%% calculate the means and the averages
freqParams.tapers = [0.5 10];
freqParams.Fs = 1000;
freqParams.dn = 0.05;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];

s_timePoints_spacing = 50;

for ei = 6:6

dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, [], [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);
end

%%
%normalization

s_length = size(s,1);
s_z = (s - repmat(s_mean,size(s,1),1))  ./ repmat(s_std,size(s,1),1);
s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.dn;


figure
imagesc(s_z')
axis xy
caxis([-3 3])


xticklabels = s_timePoints(1:s_timePoints_spacing:end);
xticks = 1:s_timePoints_spacing:length(s_timePoints);
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)

figure
plot(freq,log(mean(s)))

%need to average across all electrodes 

%%

s_timePoints_spacing = 100;
data_all_trials = concat_all_trials(Trials,experiment,driveNameAnalyze,MONKEYDIR);

for ei = 6:6

dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, [], [], [], [], freqParams.contFlag);
end

s_length = size(s,1);
s_z = (s - repmat(s_mean,size(s,1),1))  ./ repmat(s_std,size(s,1),1);
s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.dn;


figure
imagesc(s_z')
axis xy
caxis([-6 6])

xticklabels = s_timePoints(1:s_timePoints_spacing:end);
xticks = 1:s_timePoints_spacing:length(s_timePoints);
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)

%% concat all electrodes for a single trial

driveNameAnalyze = 'LM1_ECOG_3';
data_all_electrodes = concat_all_electrodes(Trials(1),experiment,driveNameAnalyze,MONKEYDIR);

[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, [], [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);

plot(freq, log(s_mean))
title('Average Spectrum Across All Electrodes for one Electrode')
xlabel('Frequency (Hz)')
ylabel('Log Amplitude')

