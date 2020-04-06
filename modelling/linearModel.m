
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
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);

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
%without zero
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
    
    trLfpData(:,trialInfo.badECoGs,:) = nan;

    
    %do target sorting
    [Pos, targids] = calcSeqReachTarg(Trials);
    Pos_Seq_1 = squeeze(Pos(:,2,:));
    trial_TargetAssigned = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);
    


    trials_index = 1:10;  
    win_sizes = [20 40 80 160 320 640];
    win_starts = -1000:100:1000;
     R_S_all = [];
    for si = 1:length(win_starts)
        for size_i = 1:length(win_sizes)
            win_start = win_starts(si);
            win_size =  win_sizes(size_i);
            
            
            feature.timePoints = (trialInfo.timeVec >= win_start) ...
                & (trialInfo.timeVec <= win_start + win_size);
            lowered_electrodes_indices = 211 + find(trialInfo.depthProfile > 0);
            %lowered_electrodes_indices = 211 + 6;
            SC32_lowered_lfp = trLfpData(:,lowered_electrodes_indices,feature.timePoints);
            X = (trLfpData(:,1:211,feature.timePoints));
            
            if win_start + win_size > max(trialInfo.timeVec)
                R_S_all(si,size_i,:) = nan(1,size(SC32_lowered_lfp,2));
                continue
            end
            
            [R_squared_prediction,A_mean] = leave_one_out_wiener(SC32_lowered_lfp(trials_index,:,:)...
                , X(trials_index,:,:));
            
            R_S_all(si,size_i,:) = mean(R_squared_prediction);
            
        end
        disp(['Finished: ' num2str(si)] )
    end
    
    figure
    imagesc(sq(mean(R_S_all,3)))
    xticklabels(win_sizes)
    xlabel('Sizes (ms)')
    yticks(1:length(win_starts))
    yticklabels(win_starts)
    title('Mean R2 for different times and window sizes')
    colorbar
    
    
    figure
    subplot(2,2,1)
    win_size = 6;
    imagesc(sq(R_S_all(:,win_size,:)))
    xlabel('Electrodes')
    ylabel('starting time')
    yticks(1:length(win_starts))
    yticklabels(win_starts)
    %caxis([0 0.8])
    colorbar
    title(['Window size ' num2str( win_sizes(win_size)) 'ms'])
    
   
    subplot(2,2,2)
    plot(win_starts,sq(R_S_all(:,win_size,:)))
    title(['Window Start Times for the same window size ' num2str( win_sizes(win_size))])
    
    
 
    subplot(2,2,3)
    win_size = 5;
    imagesc(sq(R_S_all(:,win_size,:)))
    xlabel('Electrodes')
    ylabel('starting time')
    yticks(1:length(win_starts))
    yticklabels(win_starts)
    %caxis([0 0.8])
    colorbar
    title(['Window size ' num2str( win_sizes(win_size)) 'ms'])
    
    subplot(2,2,4)
    plot(win_starts,sq(R_S_all(:,win_size,:)))
    title(['Window Start Times for the same window size ' num2str( win_sizes(win_size)) 'ms'])
    

    %neuron drop test 
    
    
end
% plot(Y(:,2))
% hold on
% plot(Y_hat(:,2))
% scatter(trialInfo.depthProfile(find(trialInfo.depthProfile > 0)), R_squared)

%%
%build lags into the system

trial_ind = 1:10;
lags = 0:1:20; %in ms
num_lags = length(lags);
%created time shifted version

lowered_electrodes_indices = 211 + find(trialInfo.depthProfile > 0);
R_squared_preds = [];
R_squared_preds_std =[];

% Use ECOG to predict SC32
for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,lowered_electrodes_indices,:);
    Y_shifted = Y(trial_ind,:,lag+1:end); 
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,1:end-lag);
    [R_squared_prediction,A_mean] = leave_one_out_wiener(Y_shifted, X_shifted);
    
    R_squared_preds(lag_i,:) = mean(R_squared_prediction);
    R_squared_preds_std(lag_i,:) = std(R_squared_prediction);
    
    lag_i/num_lags

end

figure
imagesc(R_squared_preds')
xticks(1:length(lags))
xticklabels(lags)
ylabel('SC32 Units')
xlabel('Lags (ms)')
title('Predicting SC32 using ECoG')
colorbar

%try removing the effect with common average filtering
% Use ECOG to predict SC32
R_squared_preds_car = [];
R_squared_preds_std_car =[];

for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,lowered_electrodes_indices,:);
    Y_shifted = Y(trial_ind,:,lag+1:end);
    Y_shifted_car =  Y_shifted - repmat(nanmean(Y_shifted,3),1,1,size(Y_shifted,3));
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,1:end-lag);
    
    X_shifted_car =  X_shifted - repmat(nanmean(X_shifted,3),1,1,size(X_shifted,3));
    [R_squared_prediction,A_mean] = leave_one_out_wiener(Y_shifted_car, X_shifted_car);
    
    R_squared_preds_car(lag_i,:) = mean(R_squared_prediction);
    R_squared_preds_std_car(lag_i,:) = std(R_squared_prediction);
    
    lag_i/num_lags

end

figure
imagesc(R_squared_preds_car')
xticks(1:length(lags))
xticklabels(lags)
ylabel('SC32 Units')
xlabel('Lags (ms)')
colorbar
title('CAR subtracted R2 as function of lags')


trial_ind = 1:10;
lags = 0:1:20; %in ms
num_lags = length(lags);
%created time shifted version

lowered_electrodes_indices = 211 + find(trialInfo.depthProfile > 0);
R_squared_preds_SC32 = [];
R_squared_preds_std_SC32 =[];

% Use SC32 to ECoG
for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,lowered_electrodes_indices,:);
    Y_shifted = Y(trial_ind,:,lag+1:end);
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,1:end-lag);
    [R_squared_prediction,A_mean] = leave_one_out_wiener(X_shifted, Y_shifted);
    
    R_squared_preds_SC32(lag_i,:) = mean(R_squared_prediction);
    R_squared_preds_std_SC32(lag_i,:) = std(R_squared_prediction);
    
    lag_i/num_lags

end

figure
imagesc(R_squared_preds_SC32')
xticks(1:length(lags))
xticklabels(lags)
ylabel('ECoG Units')
xlabel('Lags (ms)')
title('Predicting ECoG using Sc32')
colorbar




%try removing the effect with common average filtering
% Use ECOG to predict SC32
R_squared_preds_car_SC32 = [];
R_squared_preds_std_car_SC32 =[];

for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,lowered_electrodes_indices,:);
    Y_shifted = Y(trial_ind,:,1:end-lag);
    Y_shifted_car =  Y_shifted - repmat(nanmean(Y_shifted,3),1,1,size(Y_shifted,3));
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,lag+1:end);
    
    X_shifted_car =  X_shifted - repmat(nanmean(X_shifted,3),1,1,size(X_shifted,3));
    [R_squared_prediction,A_mean] = leave_one_out_wiener(X_shifted_car, Y_shifted_car);
    
    R_squared_preds_car_SC32(lag_i,:) = mean(R_squared_prediction);
    R_squared_preds_std_car_SC32(lag_i,:) = std(R_squared_prediction);
    
    lag_i/num_lags

end

figure
imagesc(R_squared_preds_car_SC32')
xticks(1:length(lags))
xticklabels(lags)
ylabel('ECOG Units')
xlabel('Lags (ms)')
colorbar
title('CAR subtracted R2 as a function of lags')

%use one SC32 electrodes to predict
%try removing the effect with common average filtering
% Use ECOG to predict SC32
R_squared_preds_car_SC32_e7 = [];


for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,211 + 7,:);
    Y_shifted = Y(trial_ind,:,1:end-lag);
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,lag+1:end);
    
    X_shifted_car =  X_shifted - repmat(nanmean(X_shifted,3),1,1,size(X_shifted,3));
    [R_squared_prediction,A_mean] = leave_one_out_wiener(X_shifted_car, Y_shifted);
    
    R_squared_preds_car_SC32_e7(lag_i,:) = mean(R_squared_prediction);
    
    lag_i/num_lags

end



R_squared_preds_car_SC32_e6 = [];
for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,211 + 6,:);
    Y_shifted = Y(trial_ind,:,1:end-lag);
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,lag+1:end);
    
    X_shifted_car =  X_shifted - repmat(nanmean(X_shifted,3),1,1,size(X_shifted,3));
    [R_squared_prediction,A_mean] = leave_one_out_wiener(X_shifted_car, Y_shifted);
    
    R_squared_preds_car_SC32_e6(lag_i,:) = mean(R_squared_prediction);
    
    lag_i/num_lags

end


R_squared_preds_car_SC32_e14 = [];
for lag_i = 1:num_lags

    lag = lags(lag_i);
    Y = trLfpData(:,211 + 14,:);
    Y_shifted = Y(trial_ind,:,1:end-lag);
    
    X = trLfpData(:,trialInfo.goodECoGs,:);
    X_shifted = X(trial_ind,:,lag+1:end);
    
    X_shifted_car =  X_shifted - repmat(nanmean(X_shifted,3),1,1,size(X_shifted,3));
    [R_squared_prediction,A_mean] = leave_one_out_wiener(X_shifted_car, Y_shifted);
    
    R_squared_preds_car_SC32_e14(lag_i,:) = mean(R_squared_prediction);
    
    lag_i/num_lags

end

figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
imagesc(R_squared_preds_car_SC32_e7')
xticks(1:length(lags))
xticklabels(lags)
ylabel('ECOG Units')
xlabel('Lags (ms)')
colorbar
title('SC32 E7 to model ECoG')
caxis([0 1])

subplot(1,3,2)
imagesc(R_squared_preds_car_SC32_e6')
xticks(1:length(lags))
xticklabels(lags)
ylabel('ECOG Units')
xlabel('Lags (ms)')
colorbar
title('SC32 E6 to model ECoG')
caxis([0 1])


subplot(1,3,3)
imagesc(R_squared_preds_car_SC32_e14')
xticks(1:length(lags))
xticklabels(lags)
ylabel('ECOG Units')
xlabel('Lags (ms)')
colorbar
title('SC32 E14 to model ECoG')
caxis([0 1])



