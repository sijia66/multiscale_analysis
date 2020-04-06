%%
% the goal of this
%develope the code pipeline to process lfp data
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
addpath(genpath([MONKEYDIR '/m/']))

% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180329','180329'});
Sessions = cat(1, driveSessions{:});
driveNames = Sessions(:,3);
driveNameAnalyze = {'LM1_ECOG_3'};
useSess          = ismember(driveNames, driveNameAnalyze);
DayAnalyze = Sessions(useSess,1);

% define processing information
trialInfo.sequence_random = 1;

%yes, we are filtering out the long acq time
trialInfo.filterAcq = 1;
trialInfo.filterAcq_time = 5000; %ms

trialInfo.tBefore = -1.5e3;
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;

%find the usage electrodes
trialInfo.badECoGs = [47 59 163];
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);
%trialInfo.goodSC32 = 211 + find(trialInfo.depthProfile > 0)

%define neural feature
featureROI  = [-50 200]; %in ms

%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

%% then we compare the signal obtain to that our CAR

for iD=1:nSess
    %load the trials and get the depth information
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    trialInfo.sessName =  DayAnalyze{iD};
    trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
    
    %filter out trials
    [trialInfo,Trials] = filter_trials(trialInfo,Trials);
    trialInfo.goodSC32 = 211 + find(trialInfo.depthProfile > 0)
    
    %first load an experiment trial
    expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
    load(expFN,'experiment')
    
    %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
    
    %load data
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    [trLfpData_ECOG_CAR, data_mean_ECOG] = cal_car(trLfpData(:,trialInfo.goodECoGs,:),2);
    [trLfpData_SC32_CAR,data_mean_SC32] = cal_car(trLfpData(:,trialInfo.goodSC32,:),2);

end

%%
%look at sample electrodes
ei = 1;
ti = 1;

plot(trialInfo.timeVec,sq(trLfpData(ti,ei,:)))
hold on
plot(trialInfo.timeVec,sq(data_mean_ECOG(ti,1,:)))
plot(trialInfo.timeVec,sq(trLfpData_ECOG_CAR(ti,ei,:)))
xlabel('Time from movement start (ms)')
legend({'Raw','Electrode Mean','after common average filtering'})
title(['Electrode ' num2str(ti)])

figure

ei = 24;
ti = 1;

plot(trialInfo.timeVec,sq(trLfpData(ti,ei,:)))
hold on
plot(trialInfo.timeVec,sq(data_mean_ECOG(ti,1,:)))
plot(trialInfo.timeVec,sq(trLfpData_ECOG_CAR(ti,ei,:)))
xlabel('Time from movement start (ms)')
legend({'Raw','Electrode Mean','after common average filtering'})
title(['Electrode ' num2str(ti)])

%%
% we will explore how much the activity be explained by the electrode
% variance, right?
% we assume how much variance be explained by 

window_size = 10;  %ms
num_windows = floor(length(trialInfo.timeVec) / window_size);
time_vec = (1:num_windows) * window_size;
%num_trials = size(data_mean_ECOG,1);

num_trials = 133;


var_mean = nan(num_trials,211,num_windows);

for ei = 1:211
    for ti = 1:num_trials
        
        for wi = 1:num_windows
            time_indices = ((wi-1) * window_size + 1):(wi * window_size);
            diff_vec = sq(trLfpData_ECOG_CAR(ti,ei,time_indices)) - sq(data_mean_ECOG(ti,1,time_indices));
            
            var_mean(ti,ei,wi) =  sqrt(sum(diff_vec'*diff_vec) / (length(diff_vec) - 1));
            
        end
    end
    
end
trialInfo.goodTrials = setdiff(1:num_trials,23);
plot(time_vec - trialInfo.timeReachStart,sq(nanmean(nanmean(var_mean(trialInfo.goodTrials,:,:)))))
xlabel('time from movement start(ms)')
ylabel('Standard deviation (mV)')
title({'Standard deviation over time', 'averaged all trials and electrodes'})
xline(0);


%% then look at using mean activity for prediction

win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;

targ_predict_accuracy = nan(length(win_sizes),length(start_times));
timeVec = time_vec - trialInfo.timeReachStart;

for wt = 1:length(win_sizes)
    for st_i =  1:length(start_times)
        
        timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
        
        if(start_times(st_i) + win_sizes(wt) > max(timeVec))
            continue
        end
        
        neural_feature = squeeze(mean(data_mean_ECOG(:,1,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
        targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
    end
end

imagesc(targ_predict_accuracy)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('ECOG Target prediction with mean ECoG Activity')
colorbar
