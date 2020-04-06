% need to process the trials

 trialInfo.timeReachStart = 2000;
 trialInfo.badECOG = [47,59,163];
 
depthProfile = (Trials(1).Depth{1,2})'; % in micron 2 for SC32
trialInfo.badSC32  = find(depthProfile == 0);

freqParams.tapers = [0.2 5];
freqParams.Fs = 1000;
freqParams.dn = 0.02;
freqParams.contFlag = 1;
freqParams.freqRange = [0 50];

exFeatureParams.freqROI = [0 300]; %Hz
exFeatureParams.timeROI = [1500 2500]; %ms

%%

ei = 1;
targi = 1;
tickSpacing = 10; %pixels

eval(sprintf('data = dataTarget.target%d;',targi))
eval(sprintf('dataTrials = dataTarget.target%d_trials;',targi))
dataTemp = squeeze(data(:,ei,:));
[s_z_trials,s_z_info] = caltf_trialNorm(dataTemp,dataTrials,trialInfo,freqParams);
s_timePoints = s_z_info.xlabel;
freq = s_z_info.ylabel;

time_indices = find((s_timePoints <= exFeatureParams.timeROI(2) / 1000) ...
                & (s_timePoints >= exFeatureParams.timeROI(1) / 1000)) ;
freq_indices = find(freq <= exFeatureParams.freqROI(2) ...
                & freq >= exFeatureParams.freqROI(1)) ; 

s_z_trials_mean = squeeze(mean(s_z_trials)); % average all trials

figure
imagesc(s_z_trials_mean(time_indices,freq_indices)')
colorbar
axis xy

yticklabels = freq(freq_indices);
yticklabels = round(yticklabels(1:tickSpacing:end));
yticks = linspace(1, size(s_z_trials_mean(time_indices,freq_indices), 2), ...
            numel(yticklabels));
set(gca, 'yTick', yticks, 'yTickLabel', yticklabels)
ylabel('Frequency (Hz)')

xticklabels = s_timePoints(time_indices) - 2;
xticklabels = xticklabels(1:tickSpacing:end);
xticks = linspace(1, size(s_z_trials_mean(time_indices,freq_indices), 1), ...
                numel(xticklabels));
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)
xlabel('Time from Movement Onset (s)')
caxis([-2 4])

title(sprintf('Electrode %d',ei))

%% for electrodes within a a single drive (yet to be cleaned)

% need to process the trials

timeReachStart = 2000;

freqParams.tapers = [0.2 5];
freqParams.Fs = 1000;
freqParam.dn = 0.02;        
freqParams.contFlag = 1;
freqRange = [0 300];

% extract features
preMovementTime = 400; % in ms;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_z_targ_features = [];
s_z_targ_ML = [];
s_z_targ = [];

for targi = 1:7
    
    eval(sprintf('data = dataTarget.target%d;',targi))
    eval(sprintf('dataTrials = dataTarget.target%d_trials;',targi))
    
    s_z_pool = [];
    s_z_targ_trial_features =  [];
    for ei = 1:211
        for ti = 1:size(data,1)
            
            dataTemp = squeeze(data(ti,ei,:));
            %first normalize
            dataTemp = dataTemp / norm(dataTemp);
            
            timeAqGo = dataTrials(ti).StartAq - dataTrials(ti).ReachStart + timeReachStart;
            timeGo = dataTrials(ti).Go - dataTrials(ti).ReachStart + timeReachStart; %relative to reach start
            timeTargsOn = dataTrials(ti).TargsOn - dataTrials(ti).ReachStart + timeReachStart;
            
            dataTemp = dataTemp(timeAqGo:timeTargsOn);
            
            %do frequency calculation
            dat = dataTemp';
            [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);
            freqRangeIndex = cal_index_freq(freq ,freqRange(1),freqRange(2));
            
            s_n_mean = squeeze(mean(s));
            s_n_std = squeeze(std(s));
            
            % now look at the trials
            % we are gonna get
            dataTemp = squeeze(data(ti,ei,:));
            %(timeReachStart - timeInterest  + 1) : (timeReachStart + timeInterest) ));
            %first normalize
            dataTemp = dataTemp / norm(dataTemp);
            dat = dataTemp';
            [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, [], [], [], [], freqParams.contFlag);
            
            %get z-score
            sN = size(s,1);
            s_z = (s - repmat(s_n_mean,sN,1)) ./ repmat(s_n_std,sN,1);
            
            %extract features from these plots
            s_length = size(s_z,1);
            s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.tapers(1) * 0.1;
            time_indices = find((s_timePoints <= timeReachStart / 1000) ...
                & (s_timePoints >= (timeReachStart - preMovementTime) / 1000)) ;
            freq_indices = find(freq <= freqRange(2)) ;
            s_z_preMovement = squeeze(mean(s_z(time_indices,freq_indices)));
            
            %s_z_targ_data(targi,:,:) = s_z_pool_mean;
            s_z_targ_trial_features(ti,ei,:) = s_z_preMovement;
  
        end
    end
    
    %build ML datasets
    s_z_targ_ML = cat(1,s_z_targ_ML,s_z_targ_trial_features);
    s_z_targ = [s_z_targ;ones(size(data,1),1) * targi];
    
    
    s_z_targ_features(targi,:,:) = squeeze(mean(s_z_targ_trial_features));
   
    disp(strcat('finished:',num2str(targi)))
    
end
disp('All Processing Done')


%%


% extract features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

s_z_targ_ECOG = [];  
s_z_targ_SC32 = [];

for targi = 1:7
    
    eval(sprintf('data = dataTarget.target%d;',targi))
    eval(sprintf('dataTrials = dataTarget.target%d_trials;',targi))
    
    
    %process ECOG first
    s_z_targ_trial_features = [];
    for ei = 1:size(data,2)
        if ismember(ei, trialInfo.badECOG)
            %fprintf('Ignore ECOG ET %d\n',ei)
            continue
        end
        
        dataTemp = squeeze(data(:,ei,:));
        [s_z_trials,s_z_info] = caltf_trialNorm(dataTemp,dataTrials,trialInfo,freqParams);
        s_timePoint = s_z_info.xlabel;
        freq = s_z_info.ylabel;
        
        
        %then we can extract features by averaging in time and freq domains
        % 
        time_indices = find((s_timePoints <= exFeatureParams.timeROI(2) / 1000) ...
                & (s_timePoints >= exFeatureParams.timeROI(1) / 1000)) ;
        freq_indices = find(freq <= exFeatureParams.freqROI(2) ...
                & freq >= exFeatureParams.freqROI(1)) ; 
            
        s_z_targ_trial_features(ei,:) = ...
            squeeze(mean(mean(s_z_trials(:,time_indices,freq_indices),2),3));
            
    end
    s_z_targ_ECOG = [s_z_targ_ECOG s_z_targ_trial_features] ;   
 
    %SC32
    eval(sprintf('data = dataTarget_SC32.target%d;',targi))
    eval(sprintf('dataTrials = dataTarget_SC32.target%d_trials;',targi))
    s_z_targ_trial_features = [];
    
    for ei = 1:size(data,2)
        
        if ismember(ei, trialInfo.badSC32)
            %fprintf('Ignore SC32 ET %d\n',ei)
            continue
        end
        
        dataTemp = squeeze(data(:,ei,:));
        [s_z_trials,s_z_info] = caltf_trialNorm(dataTemp,dataTrials,trialInfo,freqParams);
        s_timePoint = s_z_info.xlabel;
        freq = s_z_info.ylabel;
        
        
        %then we can extract features by averaging in time and freq domains
        % 
        time_indices = find((s_timePoints <= exFeatureParams.timeROI(2) / 1000) ...
                & (s_timePoints >= exFeatureParams.timeROI(1) / 1000)) ;
        freq_indices = find(freq <= exFeatureParams.freqROI(2) ...
                & freq >= exFeatureParams.freqROI(1)) ; 
            
        s_z_targ_trial_features(ei,:) = ...
            squeeze(mean(mean(s_z_trials(:,time_indices,freq_indices),2),3));
            
    end
    s_z_targ_SC32 = [s_z_targ_SC32 s_z_targ_trial_features] ;   
    
    fprintf('Finished %d\n',targi)
    
end

toc