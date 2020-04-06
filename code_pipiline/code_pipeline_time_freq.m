%the goal of this code is just to see if we can use LFPs directly to
%predict target direction, just like what we did in the lab


%% 
%develope the code pipeline to process lfp data
global MONKEYDIR
MONKEYDIR = 'D:/Users/Orsborn Lab/OneDrive - UW/projects/Brain EEG/data';
addpath(genpath([MONKEYDIR '/m/']))

%%
% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180408','180408'});
Sessions = cat(1, driveSessions{:});
driveNames = Sessions(:,3);
driveNameAnalyze = {'LM1_ECOG_3'};
useSess          = ismember(driveNames, driveNameAnalyze);
DayAnalyze = Sessions(useSess,1);

%%
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


freqParams.tapers = [0.2 5];
freqParams.Fs = 1000;
freqParams.dn = 0.02;
freqParams.contFlag = 1;
freqParams.freqRange = [0 200];


%define neural feature
featureROI  = [-50 200]; %in ms

%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'D:\Users\Orsborn Lab\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

%% 
win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;

s_trial_e_t_f = [];


for iD=1:nSess
    %load the trials and get the depth information
    trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
    load(trFN,'Trials')
    trialInfo.sessName =  DayAnalyze{iD};
    trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
    
    %filter out trials
    [trialInfo,Trials] = filter_trials(trialInfo,Trials);
    
    %first load an experiment trial
    expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
    load(expFN,'experiment')
    
    %load data
    [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
    
    
    for ti = 1:1
        data = squeeze(trLfpData(ti,:,:));
        [s_z_trials,s_z_info] = caltf_trial(data,freqParams);
        s_trial_e_t_f(ti,:,:,:) = s_z_trials;
        ti/size(trLfpData,1)
    end
    
    %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
    
    
    
    %calculate feature extraction
    timeVec = trialInfo.timeVec;
    
    targ_predict_accuracy = nan(length(win_sizes),length(start_times));
    targ_predict_accuracy_car = nan(length(win_sizes),length(start_times));
    
%     for wt = 1:length(win_sizes)
%         for st_i =  1:length(start_times)
%             
%             if(start_times(st_i) + win_sizes(wt) > max(timeVec))
%                 continue
%             end
%             
%             timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
%             neural_feature = squeeze(mean(trLfpData(:,1:211,timeIdx),3));
%             target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
%             targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
%             
%   
%         end
%         wt/length(win_sizes)
%     end   
end


figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,2,1)
imagesc(targ_predict_accuracy)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('Target prediction accuracy w/o CAR')
colorbar
caxis([0 1])

subplot(1,2,2)
imagesc(targ_predict_accuracy_car)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('Target prediction accuracy w/ CAR')
colorbar
caxis([0 1])

%could be due to the distribution of targets
% we also need to look at what happens to the time windows.
% can we pool together the resources. 
% we can use PCA to trace out the space 
% we should also look at the gamma band

%%
%look at the raw signal
ei = 1;
subplot(2,1,1)
for i = 1:7
    plot(trialInfo.timeVec, squeeze(mean(trLfpData(trial_TargetAssigned == i,ei,:))))
    hold on
    xlabel('Time (ms)')
    ylabel('Amp.')
    xline(0)
    xline(200)
    xline(400)
end

subplot(2,1,2)
for i = 1:7
    plot(trialInfo.timeVec, squeeze(mean(trLfpData_ECOG_CAR(trial_TargetAssigned == i,ei,:))))
    xlabel('Time (ms)')
    ylabel('Amp.')
    hold on
    xline(0)
    xline(200)
    xline(400)
end



% figure
% subplot(1,2,1)
% 
% C = confusionchart(double(trial_TargetAssigned),target_hat);
% title('Target prediction w/o ECOG CAR')
% 
% subplot(1,2,2)
% 
% C = confusionchart(double(trial_TargetAssigned),target_hat_car);
% title('Target prediction w/ ECOG CAR')
% 
% dim1 = [0.5 0.5 0.3 0.3];
% dim2 = [0.8 0.5 0.3 0.3];
% annotation('textbox',dim1, 'String',num2str(targ_predict_accuracy(1)),'FitBoxToText','on')
% annotation('textbox',dim2, 'String',num2str(targ_predict_accuracy(2)),'FitBoxToText','on')





