%the goal of this code is just to see if we can use LFPs directly to
%predict target direction, just like what we did in the lab


%% load example trials
%develope the code pipeline to process lfp data
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
%MONKEYDIR = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\data'

%load ECoG files
%selectedfile =  'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
selectedfile =  'D:\SiJia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
addpath(genpath([MONKEYDIR '/m/']))

% look at the sessions
drive_base = 'LM1_ECOG';
driveSessions = makeDriveDatabase(drive_base,{'180328','180328'});
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
%selectedfile =  'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

%% load data of a session
iD = 1;

%load the trials and get the depth information
trFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/Trials.mat'];
load(trFN,'Trials')
trialInfo.sessName =  DayAnalyze{iD};
trialInfo.depthProfile = (Trials(1).Depth{1,2})'; % in micron
trialInfo.badECoGs = [47 59 163];
trialInfo.ECoG_offset = 211;
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);
trialInfo.goodSC32 = (find( trialInfo.depthProfile > 0))';
trialInfo.goodE = [trialInfo.goodECoGs trialInfo.goodSC32+trialInfo.ECoG_offset];
trialInfo.proc_electrodes = 1:243;
trialInfo.badE = setdiff(trialInfo.proc_electrodes,trialInfo.goodE);

%filter out trials
[trialInfo,Trials] = filter_trials(trialInfo,Trials);

%first load an experiment trial
expFN = [MONKEYDIR '/' DayAnalyze{iD} '/' Trials(1).Rec '/rec' Trials(1).Rec '.experiment.mat'];
load(expFN,'experiment')

%get the chamber information
P_chamber_ECOG = drmap_coordTransform_microdrive2chamber(experiment,1);
P_chamber_ECOG = P_chamber_ECOG';

P_chamber_SC32 = drmap_coordTransform_microdrive2chamber(experiment,2);
P_chamber_SC32 = P_chamber_SC32';

%only look at good electrodes
P_chamber_ECoG = P_chamber_ECOG(trialInfo.goodECoGs,:);
P_chamber_SC32 = P_chamber_SC32(trialInfo.goodSC32,:);

%get car filter info
%get common average referencing filters
refMatrix_car_ECoG = ...
    car_matrix(P_chamber_ECoG);

refMatrix_car_SC32 = ...
    car_matrix(P_chamber_SC32);

%get laplacian matrix
distThres = 0.75;
[refMatrix_lap_ECoG, P_chamber_lap_ECoG,electrode_index_ECoG] = ...
                laplacian_matrix(P_chamber_ECoG,distThres);
distThres = 1.5;
[refMatrix_lap_SC32, P_chamber_lap_SC32,electrode_index_SC32] = ...
                laplacian_matrix(P_chamber_SC32,distThres);

%load data
[trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR);
trLfpData_ECOG = trLfpData(:,trialInfo.goodECoGs,:);
trLfpData_SC32 = trLfpData(:,trialInfo.goodSC32,:);

X = tensor(trLfpData(:,trialInfo.goodECoGs,:));
trLfpData_ECOG_CAR = double(ttm(X,refMatrix_car_ECoG,2));
trLfpData_ECOG_lap = double(ttm(X,refMatrix_lap_ECoG,2));

X = tensor(trLfpData(:,trialInfo.goodSC32,:));
trLfpData_SC32_CAR = double(ttm(X,refMatrix_car_SC32,2));
trLfpData_SC32_lap = double(ttm(X,refMatrix_lap_SC32,2));
%trLfpData_ECOG_CAR = cal_car(trLfpData(:,trialInfo.ECOG_indices,:),2);
%trLfpData_SC32_CAR = cal_car(trLfpData(:,211 + find(trialInfo.depthProfile > 0),:),2);

%do target sorting
[Pos, targids] = calcSeqReachTarg(Trials);
Pos_Seq_1 = squeeze(Pos(:,2,:));
trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
%% do LDA for 

%% do LDA

win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;
%calculate feature extraction
timeVec = trialInfo.timeVec;

targ_predict_accuracy = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_SC32 = nan(length(win_sizes),length(start_times));

max_accuracy_SC32 = 0;
max_accuracy_ECOG = 0;
disp('only use length(SC32) - 1 SC32 electrodes')

for wt = 1:length(win_sizes)
    for st_i =  1:length(start_times)
        
        timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
        if(start_times(st_i) + win_sizes(wt) > max(timeVec))
            continue
        end
        
        neural_feature = squeeze(mean(trLfpData_ECOG(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
        targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
        
        if  targ_predict_accuracy(wt,st_i) > max_accuracy_ECOG
            target_hat_max_ECOG =  target_hat;
            max_accuracy_ECOG = targ_predict_accuracy(wt,st_i);
        end
        
    
        
        %SC32
        neural_feature = squeeze(mean(trLfpData_SC32(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
        targ_predict_accuracy_SC32(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
        if  targ_predict_accuracy_SC32(wt,st_i) > max_accuracy_SC32
            target_hat_max_SC32 =  target_hat;
            max_accuracy_SC32 = targ_predict_accuracy_SC32(wt,st_i);
        end
       
        
    end
    wt/length(win_sizes)
end


%% do LDA

win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;
%calculate feature extraction
timeVec = trialInfo.timeVec;

targ_predict_accuracy = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_car = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_SC32 = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_SC32_car = nan(length(win_sizes),length(start_times));

max_accuracy_SC32 = 0;
max_accuracy_ECOG = 0;
disp('only use length(SC32) - 1 SC32 electrodes')

for wt = 1:length(win_sizes)
    for st_i =  1:length(start_times)
        
        timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
        if(start_times(st_i) + win_sizes(wt) > max(timeVec))
            continue
        end
        
        neural_feature = squeeze(mean(trLfpData_ECOG(:,:,timeIdx),3));
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
        neural_feature = squeeze(mean(trLfpData_SC32(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
        targ_predict_accuracy_SC32(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
        if  targ_predict_accuracy_SC32(wt,st_i) > max_accuracy_SC32
            target_hat_max_SC32 =  target_hat;
            max_accuracy_SC32 = targ_predict_accuracy_SC32(wt,st_i);
        end
        
        %SC32 CAR
        %try
            neural_feature = squeeze(mean(trLfpData_SC32_CAR(:,1:(length(trialInfo.goodSC32) - 1),timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            targ_predict_accuracy_SC32_car(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        %catch
            %continue
            
        %end
        neural_feature = squeeze(mean(trLfpData_SC32_lap(:,1:(length(trialInfo.goodSC32) - 2),timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
        targ_predict_accuracy_SC32_lap(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
    end
    wt/length(win_sizes)
end
%% visualize the results
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
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

subplot(1,3,2)
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

subplot(1,3,3)
imagesc(targ_predict_accuracy_lap)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('ECOG Target prediction accuracy w/ Laplacian')
colorbar
caxis([0 1])
saveas(gca, 'compare_rereferencing_tech_ECoG.png')


figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
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

subplot(1,3,2)
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

subplot(1,3,3)
imagesc(targ_predict_accuracy_SC32_lap)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('SC32 Target prediction accuracy w/ Laplacian')
colorbar
caxis([0 1])

saveas(gca, 'compare_rereferencing_tech_SC32.png')
%% neuron drop test
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
%% neuron drop test for SC32

%choose one timewindow and start time
st_i = 7;
wt = 8;
start_time = 200; %ms
win_size = 400; %ms
timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));

num_units = length(trialInfo.goodSC32);
targ_predict_accuracy_SC32_neuron = nan(num_units,1);

neural_feature = squeeze(mean(trLfpData_SC32(:,:,timeIdx),3));
target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
full_accuracy_SC32_neuron ...
    = sum(target_hat == trial_TargetAssigned) / length(target_hat);


%SC32
for u_i = 1:num_units
    
    use_for_training = true(num_units,1);
    %drop the u_ith neuron
    use_for_training(u_i) = false;
    
    %use the rest for training
    neural_feature = squeeze(mean(trLfpData_SC32(:,use_for_training,timeIdx),3));
    try
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
    catch
        continue
    end
    targ_predict_accuracy_SC32_neuron(u_i) ...
        = sum(target_hat == trial_TargetAssigned) / length(target_hat);
    
end

figure
scatter(trialInfo.goodSC32, targ_predict_accuracy_SC32_neuron)
ylabel('Accuracy')
xlabel('Dropped units')
ylim([0 0.5])
saveas(gca, 'neuron drop test.png')

fig = figure;
[map,mapInfo]= scatterToMatrix(P_chamber_SC32,targ_predict_accuracy_SC32_neuron);
imagesc(map)
scatterToMatrix_plot_label(fig,map,mapInfo)
%% examine the result of neuron drop test
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
%% look at the raw signal
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
%% look at confusion mat 
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
%% look at frequency response separately for ECoG
load('filter_bank.mat')
%now we should have the frequency bands and so on.

win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;
%calculate feature extraction
timeVec = trialInfo.timeVec;

targ_predict_accuracy = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_frequency = nan(length(filter_bank), length(win_sizes),length(start_times));
targ_predict_accuracy_SC32_frequency = nan(length(filter_bank), length(win_sizes),length(start_times));

max_accuracy_SC32 = 0;
max_accuracy_ECOG = 0;
disp('only use length(SC32) - 1 SC32 electrodes')

for wt = 1:length(win_sizes)
    for st_i =  1:length(start_times)
        
        timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
        if(start_times(st_i) + win_sizes(wt) > max(timeVec))
            continue
        end
        
        neural_feature = squeeze(mean(trLfpData_ECOG(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
        targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
        
        if  targ_predict_accuracy(wt,st_i) > max_accuracy_ECOG
            target_hat_max_ECOG =  target_hat;
            max_accuracy_ECOG = targ_predict_accuracy(wt,st_i);
        end
    end
    disp(wt / length(win_sizes))
end

for fi = 1:length(filter_bank)
    %exchange the trial and time dimensions
    trLfpData_temp  = permute(trLfpData_ECOG, [3 2 1]);
    trLfpData_temp  = filter(filter_bank(fi),trLfpData_temp);  %1 for the time dimension
    trLfpData_temp  = permute(trLfpData_temp, [3 2 1]);
    
    for wt = 1:length(win_sizes)
        for st_i =  1:length(start_times)
            %calculate time indix to average
            timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
            if(start_times(st_i) + win_sizes(wt) > max(timeVec))
                continue
            end
            
            neural_feature = squeeze(mean(trLfpData_temp(:,:,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
            targ_predict_accuracy_frequency(fi,wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        end
    end
    disp(fi / length(filter_bank))
end

%now do SC32
disp('now processing SC32')
for fi = 1:length(filter_bank)
    %exchange the trial and time dimensions
    trLfpData_temp  = permute(trLfpData_SC32, [3 2 1]);
    trLfpData_temp  = filter(filter_bank(fi),trLfpData_temp);  %1 for the time dimension
    trLfpData_temp  = permute(trLfpData_temp, [3 2 1]);
    
    for wt = 1:length(win_sizes)
        for st_i =  1:length(start_times)
            %calculate time indix to average
            timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
            if(start_times(st_i) + win_sizes(wt) > max(timeVec))
                continue
            end
            
            neural_feature = squeeze(mean(trLfpData_temp(:,:,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            targ_predict_accuracy_SC32_frequency(fi,wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        end
    end
    disp(fi / length(filter_bank))
end

save('target_predict_accuracies',...
     'targ_predict_accuracy',...
     'targ_predict_accuracy_SC32',...
     'targ_predict_accuracy_frequency',...
     'targ_predict_accuracy_SC32_frequency'...
     );
%% plot out examples
plot_step_step = 2;

figure('units','normalized','outerposition',[0 0 0.6 0.5])
subplot(1,2,1)
imagesc(targ_predict_accuracy)
ylabel('Window sizes (s)')
xlabel('Start times (s)')
xticks(1:plot_step_step:length(start_times))
xticklabels(start_times(1:plot_step_step:end) / 1000)
yticks(1:plot_step_step:length(win_sizes))
yticklabels(win_sizes(1:plot_step_step:end)/1000)
title('ECOG broadband')
colorbar
caxis([0 0.8])
set(gca, 'fontSize', 14)

subplot(1,2,2)
imagesc(targ_predict_accuracy_SC32)
ylabel('Window sizes (s)')
xlabel('Start times (s)')
xticks(1:plot_step_step:length(start_times))
xticklabels(start_times(1:plot_step_step:end)/1000)
yticks(1:plot_step_step:length(win_sizes))
yticklabels(win_sizes(1:plot_step_step:end)/1000)
title('SC32 Broadband')
colorbar
caxis([0 0.8])
set(gca, 'fontSize', 14)

saveas(gca, 'ECoG_SC32_example_decoding accuracy.png')

%% get the max
targ_predict_accuracy_frequency_max(1) = max(targ_predict_accuracy,[],'all');
targ_predict_accuracy_frequency_max(2:5) = ...
    max(targ_predict_accuracy_frequency,[],[2 3]);

targ_predict_accuracy_SC32_frequency_max(1) = max(targ_predict_accuracy_SC32,[],'all');
targ_predict_accuracy_SC32_frequency_max(2:5)  = ...
    max(targ_predict_accuracy_SC32_frequency,[],[2 3]);

figure
plot(targ_predict_accuracy_frequency_max, 'b','LineWidth',2)
hold on
plot(targ_predict_accuracy_SC32_frequency_max, 'y','LineWidth',2)
ylabel('Accuracy of target decoding')

set(gca,'xtick', 1:5)
set(gca,'xticklabel',...
    {'Broadband',...
        [num2str(freq_bands(1,1)),'-',num2str(freq_bands(1,2)),' Hz'],...
        [num2str(freq_bands(2,1)),'-',num2str(freq_bands(2,2)),' Hz'],...
        [num2str(freq_bands(3,1)),'-',num2str(freq_bands(3,2)),' Hz'],...
        [num2str(freq_bands(4,1)),'-',num2str(freq_bands(4,2)),' Hz'],...
        })
legend('ECoG', 'LFP')
xtickangle(-30)
set(gca, 'fontSize', 18)
saveas(gca, 'max decoding accuracy.png')

 %% plot out the results ECoG
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,5,1)
imagesc(targ_predict_accuracy)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('ECOG unfiltered')
colorbar
caxis([0 1])

for fi = 1:length(filter_bank)
    subplot(1,5,fi+1)
    imagesc(squeeze(targ_predict_accuracy_frequency(fi,:,:)))
    ylabel('Window sizes (ms)')
    xlabel('Start times (ms)')
    xticks(1:length(start_times))
    xticklabels(start_times)
    yticks(1:length(win_sizes)) 
    yticklabels(win_sizes)
    title(['ECoG:',num2str(freq_bands(fi,1)),' - ',num2str(freq_bands(fi,2)), ' Hz'])
    colorbar
    caxis([0 1])    
end
saveas(gca, 'ECoG decoding accuracy.png')
%% SC32 across frequencies

figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,5,1)
imagesc(targ_predict_accuracy_SC32)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('SC32 broadband')
colorbar
caxis([0 1])

for fi = 1:length(filter_bank)
    subplot(1,5,fi+1)
    imagesc(squeeze(targ_predict_accuracy_SC32_frequency(fi,:,:)))
    ylabel('Window sizes (ms)')
    xlabel('Start times (ms)')
    xticks(1:length(start_times))
    xticklabels(start_times)
    yticks(1:length(win_sizes))
    yticklabels(win_sizes)
    title(['SC32:',num2str(freq_bands(fi,1)),' - ',num2str(freq_bands(fi,2)), ' Hz'])
    colorbar
    caxis([0 1])    
end

saveas(gca, 'SC32 decoding accuracy.png')
%% look at LDA for filtered signals with rms 

load('filter_bank.mat')
%now we should have the frequency bands and so on.

win_sizes = [10 20 40 80 100 160 200 400 500 750 1000 1250 1500 1750 2000];
start_times = -1000:200:1000;
%calculate feature extraction
timeVec = trialInfo.timeVec;

targ_predict_accuracy = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_SC32 = nan(length(win_sizes),length(start_times));
targ_predict_accuracy_frequency = nan(length(filter_bank), length(win_sizes),length(start_times));
targ_predict_accuracy_SC32_frequency = nan(length(filter_bank), length(win_sizes),length(start_times));

max_accuracy_SC32 = 0;
max_accuracy_ECOG = 0;
disp('only use length(SC32) - 1 SC32 electrodes')

for wt = 1:length(win_sizes)
    for st_i =  1:length(start_times)
        
        timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
        if(start_times(st_i) + win_sizes(wt) > max(timeVec))
            continue
        end
        
        neural_feature = squeeze(rms(trLfpData_ECOG(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
        targ_predict_accuracy(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
        
        if  targ_predict_accuracy(wt,st_i) > max_accuracy_ECOG
            target_hat_max_ECOG =  target_hat;
            max_accuracy_ECOG = targ_predict_accuracy(wt,st_i);
        end
        
        neural_feature = squeeze(rms(trLfpData_SC32(:,:,timeIdx),3));
        target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
        targ_predict_accuracy_SC32(wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        
        
    end
    disp(wt / length(win_sizes))
end

for fi = 1:length(filter_bank)
    %exchange the trial and time dimensions
    trLfpData_temp  = permute(trLfpData_ECOG, [3 2 1]);
    trLfpData_temp  = filter(filter_bank(fi),trLfpData_temp);  %1 for the time dimension
    trLfpData_temp  = permute(trLfpData_temp, [3 2 1]);
    
    for wt = 1:length(win_sizes)
        for st_i =  1:length(start_times)
            %calculate time indix to average
            timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
            if(start_times(st_i) + win_sizes(wt) > max(timeVec))
                continue
            end
            
            neural_feature = squeeze(rms(trLfpData_temp(:,:,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature(:,1:64),trial_TargetAssigned);
            targ_predict_accuracy_frequency(fi,wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        end
    end
    disp(fi / length(filter_bank))
end

%now do SC32
disp('now processing SC32')
for fi = 1:length(filter_bank)
    %exchange the trial and time dimensions
    trLfpData_temp  = permute(trLfpData_SC32, [3 2 1]);
    trLfpData_temp  = filter(filter_bank(fi),trLfpData_temp);  %1 for the time dimension
    trLfpData_temp  = permute(trLfpData_temp, [3 2 1]);
    
    for wt = 1:length(win_sizes)
        for st_i =  1:length(start_times)
            %calculate time indix to average
            timeIdx = (timeVec >= start_times(st_i)) & (timeVec <=  start_times(st_i) + win_sizes(wt));
            if(start_times(st_i) + win_sizes(wt) > max(timeVec))
                continue
            end
            
            neural_feature = squeeze(rms(trLfpData_temp(:,:,timeIdx),3));
            target_hat = runLeaveOneOutClassification(neural_feature,trial_TargetAssigned);
            targ_predict_accuracy_SC32_frequency(fi,wt,st_i) = sum(target_hat == trial_TargetAssigned) / length(target_hat);
        end
    end
    disp(fi / length(filter_bank))
end

save('target_predict_accuracies_rms',...
     'targ_predict_accuracy',...
     'targ_predict_accuracy_SC32',...
     'targ_predict_accuracy_frequency',...
     'targ_predict_accuracy_SC32_frequency'...
     );

%% plot out the results ECoG using rms
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,5,1)
imagesc(targ_predict_accuracy)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('ECOG unfiltered')
colorbar
caxis([0 1])

for fi = 1:length(filter_bank)
    subplot(1,5,fi+1)
    imagesc(squeeze(targ_predict_accuracy_frequency(fi,:,:)))
    ylabel('Window sizes (ms)')
    xlabel('Start times (ms)')
    xticks(1:length(start_times))
    xticklabels(start_times)
    yticks(1:length(win_sizes))
    yticklabels(win_sizes)
    title(['ECoG:',num2str(freq_bands(fi,1)),' - ',num2str(freq_bands(fi,2)), ' Hz'])
    colorbar
    caxis([0 1])    
end
saveas(gca, 'ECoG decoding accuracy_rms.png')

%% plot out SC32 
figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,5,1)
imagesc(targ_predict_accuracy_SC32)
ylabel('Window sizes (ms)')
xlabel('Start times (ms)')
xticks(1:length(start_times))
xticklabels(start_times)
yticks(1:length(win_sizes))
yticklabels(win_sizes)
title('ECOG unfiltered')
colorbar
caxis([0 1])

for fi = 1:length(filter_bank)
    subplot(1,5,fi+1)
    imagesc(squeeze(targ_predict_accuracy_SC32_frequency(fi,:,:)))
    ylabel('Window sizes (ms)')
    xlabel('Start times (ms)')
    xticks(1:length(start_times))
    xticklabels(start_times)
    yticks(1:length(win_sizes))
    yticklabels(win_sizes)
    title(['SC32:',num2str(freq_bands(fi,1)),' - ',num2str(freq_bands(fi,2)), ' Hz'])
    colorbar
    caxis([0 1])    
end

saveas(gca, 'SC32 decoding accuracy_rms.png')


%% equal number of contacts LDA
win_sizes = [100 200 400];
start_times = -1000:200:400;
%calculate feature extraction
timeVec = trialInfo.timeVec;

trLfpData_sample = trLfpData_ECOG(:,1:64,:);
targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
                                        win_sizes,start_times);
max(targ_predict_accuracy,[], 'all')

%% sequence scan the space
%to scan the electrode space
NUM_E_SAMPL = 14; %same as the LFP electrodes
num_intervals = floor(length(trialInfo.goodECoGs) / NUM_E_SAMPL);
%essentially doing a matricization of the recordings
test_intervals = reshape(trialInfo.goodECoGs(1:NUM_E_SAMPL*num_intervals),...
                                             num_intervals, NUM_E_SAMPL);                                      
%now num_intervals by NUM_E_SAMPL                                    
test_intervals = test_intervals';
targ_acc_intervals = nan(num_intervals,1);

%iterate over the samples
for i = 1:num_intervals
    trLfpData_sample = trLfpData_ECOG(:,test_intervals(i,:),:);
    targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
                                        win_sizes,start_times);
    targ_acc_intervals(i) = max(targ_predict_accuracy,[], 'all');
end

disp('Finished interval LDA')
figure 
scatter(1:num_intervals, targ_acc_intervals)

%% figure out the centroid of each cluster
%figure out the averaging matrix 
matrix_avg = zeros(num_intervals, length(trialInfo.goodECoGs));

for i = 1:num_intervals
    matrix_avg(i,test_intervals(i,:)) = 1/NUM_E_SAMPL;
end
P_chamber_clusters =  matrix_avg * P_chamber_ECoG;


%% divide the regions into num_intervals mutually exclusive regions
[idx,C] = kmeans(P_chamber_ECoG(1:NUM_E_SAMPL*num_intervals,:),num_intervals);
%build the cluster matrix
rng('default')  % For reproducibility
cluster_avg = zeros(num_intervals, length(P_chamber_ECoG));
for i = 1:num_intervals
    cluster_avg(i,find(idx == i)) = 1/NUM_E_SAMPL;
end
P_chamber_clusters =  cluster_avg * P_chamber_ECoG;
scatter(P_chamber_clusters(:,1),P_chamber_clusters(:,2), ...
    [],'filled');
figure
gscatter(P_chamber_ECoG(1:NUM_E_SAMPL*num_intervals,1),...
         P_chamber_ECoG(1:NUM_E_SAMPL*num_intervals,2),...
         idx)
%% 
 % write our own partitioning algoriths by greedy seearch
 NUM_E_SAMPL;
 P_chamber = P_chamber_ECoG;
 P_chamber_part = nan(length(P_chamber),1);
 unclustered_ind =  1:length(P_chamber);
 cluster_ind = nan(num_intervals,NUM_E_SAMPL);
 clustered = [];
 
 for i = 1:num_intervals
     %pick the first recording site
     p_ind =  unclustered_ind(1);
     
     n_clustered = size(unclustered_ind,1);
     
     vec_diff = P_chamber(unclustered_ind,:)...
                           -repmat(P_chamber(p_ind,:),n_clustered,1);
     
     dist_p_rest = vecnorm(vec_diff,2,2);
     [B,I] = sort(dist_p_rest);
     
     cluster_ind(i,:) = unclustered_ind(I(1:NUM_E_SAMPL));
     P_chamber_part(unclustered_ind(I(1:NUM_E_SAMPL))) = i;
     
     %prepare for next iteration, delete the clustered points
     unclustered_ind(I(1:NUM_E_SAMPL)) = [];

                       
 end
 %build the mean matrix.
 matrix_avg = zeros(num_intervals, length(trialInfo.goodECoGs));

for i = 1:num_intervals
    matrix_avg(i,cluster_ind(i,:)) = 1/NUM_E_SAMPL;
end
P_chamber_clusters =  matrix_avg * P_chamber_ECoG;
 
 gscatter(P_chamber_ECoG(:,1),...
         P_chamber_ECoG(:,2),...
         P_chamber_part)
%% use ths clustering algorithm to cal lDA
%calculate the SC32 accuracy

trLfpData_sample = trLfpData_SC32;
targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
    win_sizes,start_times);
targ_acc_intervals_SC32 = max(targ_predict_accuracy,[], 'all');

NUM_ECoG = 64;
trLfpData_sample = trLfpData_ECOG;
targ_predict_accuracy = lda_sample_time(trLfpData_sample(:,1:NUM_ECoG,:),timeVec,trial_TargetAssigned, ...
    win_sizes,start_times);
targ_acc_intervals_ECoG = max(targ_predict_accuracy,[], 'all');

NUM_ECoG = 107;
trLfpData_sample = trLfpData_ECOG;
targ_predict_accuracy = lda_sample_time(trLfpData_sample(:,1:NUM_ECoG,:),timeVec,trial_TargetAssigned, ...
    win_sizes,start_times);
targ_acc_intervals_ECoG_max = max(targ_predict_accuracy,[], 'all');

targ_acc_intervals = nan(num_intervals,1);
%iterate over the samples
for i = 1:num_intervals
    trLfpData_sample = trLfpData_ECOG(:,cluster_ind(i,:),:);
    targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
                                        win_sizes,start_times);
    targ_acc_intervals(i) = max(targ_predict_accuracy,[], 'all');
end

disp('Finished interval LDA')

%report a a histogram
figure
histogram(targ_acc_intervals)
hold on
xline(targ_acc_intervals_SC32, '','SC32 14 electrodes','LineWidth',2)
xline(targ_acc_intervals_ECoG, '','ECoG 64 electrodes','LineWidth',2)
xline(targ_acc_intervals_ECoG_max, '','ECoG 107 electrodes','LineWidth',2)
ylabel('Counts')
xlabel('Accuracy in decoding target directions')
set(gca,'fontSize',18)
saveas(gca,'compare lda same dimension_histogram.png')
%% plot comparision in the spatial grid
full_circle_size = 1000;
bubble_sizes = (targ_acc_intervals  - min(targ_acc_intervals)) ...
                / (max(targ_acc_intervals) - (min(targ_acc_intervals))) + 0.001;
bubble_size_SC32 = (targ_acc_intervals_SC32 - min(targ_acc_intervals)) ...
                / (max(targ_acc_intervals) - (min(targ_acc_intervals))) + 0.001;
figure
data1 = gscatter(P_chamber_ECoG(:,1),...
         P_chamber_ECoG(:,2),...
         P_chamber_part)
hold on
    scatter(P_chamber_clusters(:,1),P_chamber_clusters(:,2), ...
        full_circle_size*bubble_sizes,'r','filled');

    scatter(mean(P_chamber_SC32(:,1)),mean(P_chamber_SC32(:,2)), ...
        full_circle_size*bubble_size_SC32,'b','filled');
lgnd = legend([data1], 'Location', 'eastoutside');
xlabel('mm')
ylabel('mm')
saveas(gca,'compare lda same dimension_scatter.png')

%% calculate the accuracy for the hub electrodes
% neuron drop test for SC32 without the hub electrodes
win_sizes = [100 200 400];
start_times = [0 200];
%calculate feature extraction
trialInfo.hubEs = [7 17 20];
trialInfo.otherEs  = setdiff(trialInfo.goodSC32,trialInfo.hubEs);

timeVec = trialInfo.timeVec;

trLfpData_sample = trLfpData(:,trialInfo.hubEs + trialInfo.ECoG_offset,:);
%trLfpData_sample = trLfpData_SC32;
targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
                                        win_sizes,start_times);
hub_accuracy = max(targ_predict_accuracy,[], 'all');

trLfpData_sample = trLfpData_SC32;
targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
                                        win_sizes,start_times);
SC32_accuracy = max(targ_predict_accuracy,[], 'all');

num_iterations =  10;
test_num_ECOG = 1:length(trialInfo.otherEs);
max_num_units = length(trialInfo.otherEs);

targ_predict_accuracy_drop = nan(length(test_num_ECOG),num_iterations);

rng('default')  % For reproducibility

for iter_i = 1:num_iterations
    for u_i = 1:length(test_num_ECOG)
        
        rand_trial_order = randperm(max_num_units);
        use_for_training = rand_trial_order(1:test_num_ECOG(u_i));
        %get the actual electrode number
        use_for_training = trialInfo.otherEs(use_for_training);
        
        trLfpData_sample = trLfpData(:,use_for_training + trialInfo.ECoG_offset,:);
        %trLfpData_sample = trLfpData_SC32;
        targ_predict_accuracy = lda_sample_time(trLfpData_sample,timeVec,trial_TargetAssigned, ...
                                        win_sizes,start_times);
                                    
                                    
        targ_predict_accuracy_drop(u_i,iter_i) = max(targ_predict_accuracy, [], 'all');
    end
    iter_i / num_iterations
end

figure
boxplot(targ_predict_accuracy_drop')
xlabel('Number of non-hub intracortical electrodes')
ylabel('Target classification accuracy')
yline(hub_accuracy,'--r' ,'LineWidth',2)
yline(SC32_accuracy,'--r' ,'LineWidth',2)
set(gca, 'fontSize', 18)
saveas(gca, 'SC32 neuron drop test.png')



