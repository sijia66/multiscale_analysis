%the goal of this code is just to see if we can use LFPs directly to
%predict target direction, just like what we did in the lab


%% 
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
trialInfo.badECoGs = [47 59 163];
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

%% 
win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;

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
    trLfpData_ECOG_CAR = cal_car(trLfpData(:,trialInfo.ECOG_indices,:),2);
    trLfpData_SC32_CAR = cal_car(trLfpData(:,211 + find(trialInfo.depthProfile > 0),:),2);
    
    %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
    
    
    
    %calculate feature extraction
    timeVec = trialInfo.timeVec;
    
    targ_predict_accuracy = nan(length(win_sizes),length(start_times));
    targ_predict_accuracy_car = nan(length(win_sizes),length(start_times));
    
   max_accuracy_SC32 = 0;
   max_accuracy_ECOG = 0;
    
    for wt = 1:length(win_sizes)
        for st_i =  1:length(start_times)
            
            timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
            if(start_times(st_i) + win_sizes(wt) > max(timeVec))
                continue
            end
            

            neural_feature = squeeze(mean(trLfpData(:,1:64,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
            targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
            
            if  targ_predict_accuracy(wt,st_i) > max_accuracy_ECOG
               target_hat_max_ECOG =  target_hat;
               max_accuracy_ECOG = targ_predict_accuracy(wt,st_i);
            end
            
            %CAR
            neural_feature = squeeze(mean(trLfpData_ECOG_CAR(:,:,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
            targ_predict_accuracy_car(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
            
            %SC32
            neural_feature = squeeze(mean(trLfpData(:,211 + find(trialInfo.depthProfile > 0),timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            targ_predict_accuracy_SC32(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
            
            if  targ_predict_accuracy_SC32(wt,st_i) > max_accuracy_SC32
                 target_hat_max_SC32 =  target_hat;
                 max_accuracy_SC32 = targ_predict_accuracy_SC32(wt,st_i);
            end
            
            %SC32 CAR
            neural_feature = squeeze(mean(trLfpData_SC32_CAR(:,:,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            targ_predict_accuracy_SC32_car(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
  
        end
        wt/length(win_sizes)
    end
    
    %neuron drop test
    %choose one timewindow and start time
    start_time = 200; %ms
    win_size = 400; %ms
    timeIdx = (timeVec >= start_time & (timeVec <=  start_time + win_size));
    max_num_units = 211;
    num_units = 211;
    unit_step = 4;
    
    num_iterations = 10;
    test_num_ECOG = 1:unit_step:num_units;
    targ_predict_accuracy_ECOG_neuron = nan(num_iterations,length(test_num_ECOG));
    %ECOG
    for iter_i = 1:num_iterations
        for u_i = 1:length(test_num_ECOG)
            
            use_for_training = false(max_num_units,1);
            rand_trial_order = randperm(max_num_units);
            use_for_training(rand_trial_order(1:test_num_ECOG(u_i))) = true;
            
            neural_feature = squeeze(mean(trLfpData(:,use_for_training,timeIdx),3));
            
            try
                target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            catch
                continue
            end
            targ_predict_accuracy_ECOG_neuron(iter_i,u_i) ...
                = sum(target_hat == trial_TargetAssigned) / length(target_hat);
            
            
        end
    end
    
    %repeat this for SC32
    max_num_units = sum(trialInfo.depthProfile > 0);
    num_units = max_num_units;
    unit_step = 1;

    test_num_SC32 = 1:unit_step:num_units;
    targ_predict_accuracy_SC32_neuron = nan(num_iterations,length(test_num_SC32));

    %start the same sequence
    rng default
    
    for iter_i = 1:num_iterations
        for u_i = 1:length(test_num_SC32)
            
            use_for_training = false(max_num_units,1);
            rand_trial_order = randperm(max_num_units);
            use_for_training(rand_trial_order(1:test_num_SC32(u_i))) = true;
            
            neural_feature = squeeze(mean(trLfpData(:,use_for_training,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            targ_predict_accuracy_SC32_neuron(iter_i,u_i) ...
                = sum(target_hat == trial_TargetAssigned) / length(target_hat);
            
            
        end
    end
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
title('ECOG Target prediction accuracy w/o CAR')
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
title('ECOG Target prediction accuracy w/ CAR')
colorbar
caxis([0 1])



figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,2,1)
imagesc(targ_predict_accuracy_SC32)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('SC32 Target prediction accuracy w/o CAR')
colorbar
caxis([0 1])

subplot(1,2,2)
imagesc(targ_predict_accuracy_SC32_car)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('SC32 Target prediction accuracy w/ CAR')
colorbar
caxis([0 1])

%examine the result of neuron drop test
figure
plot(test_num_ECOG,targ_predict_accuracy_ECOG_neuron)
xlabel('Number of Units (randomized)')
ylabel('Validated Accuracy')
title('Unit drop test for ECOG')

figure
plot(test_num_SC32,targ_predict_accuracy_SC32_neuron)
xlabel('Number of Units (randomized)')
ylabel('Validated Accuracy')
title('Unit drop test for SC32')



%could be due to the distribution of targets
% we also need to look at what happens to the time windows.
% we can use PCA to trace out the space 
% we should also look at the gamma band

%%
%look at the raw signal
ei = 1;
figure
subplot(2,1,1)
for i = 1:7
    plot(trialInfo.timeVec, squeeze(mean(trLfpData(trial_TargetAssigned == i,ei,:))))
    hold on
    xlabel('Time (ms)')
    ylabel('Amp.')
    title(['ECOG Example Electrode ' num2str(ei)])
    xline(0)
    xline(200)
    xline(400)
end

subplot(2,1,2)
for i = 1:7
    plot(trialInfo.timeVec, squeeze(mean(trLfpData_ECOG_CAR(trial_TargetAssigned == i,ei,:))))
    xlabel('Time (ms)')
    ylabel('Amp.')
    title(['ECOG Example Electrode ' num2str(ei) ' CAR'])
    hold on
    xline(0)
    xline(200)
    xline(400)
end


ei = 6;
figure
subplot(2,1,1)
for i = 1:7
    plot(trialInfo.timeVec, squeeze(mean(trLfpData(trial_TargetAssigned == i,ei + 211,:))))
    hold on
    xlabel('Time (ms)')
    ylabel('Amp.')
    title(['SC32 Example Electrode ' num2str(ei)])
    xline(0)
    xline(200)
    xline(400)
end

subplot(2,1,2)
for i = 1:7
    plot(trialInfo.timeVec, squeeze(mean(trLfpData_SC32_CAR(trial_TargetAssigned == i,1,:)))) %arbitrary
    xlabel('Time (ms)')
    ylabel('Amp.')
    title(['SC32 Example Electrode ' num2str(ei) ' CAR'])
    hold on
    xline(0)
    xline(200)
    xline(400)
end


%%
%look at confusion mat 
figure
subplot(1,2,1)
cm = confusionchart(double(trial_TargetAssigned),target_hat_max_ECOG)
title('ECOG')
set(gca, 'fontSize',28)

subplot(1,2,2)
cm = confusionchart(double(trial_TargetAssigned),target_hat_max_SC32)
title('SC32')
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
set(gca, 'fontSize',28)




