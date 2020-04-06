%search for the time period that has greater mutual information with time

timeReachStart = 2000;
timeInterest = 300;

freqParams.tapers = [0.2 5];
freqParams.Fs = 1000;
freqParam.dn = 0.02;
freqParams.contFlag = 1;
freqRange = [0 300];

% extract features
ei = 178;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




s_z_targ_features = [];
s_z_targ_ML = [];
s_z_targ = [];

for targi = 1:7
    
    eval(sprintf('data = dataTarget.target%d;',targi))
    eval(sprintf('dataTrials = dataTarget.target%d_trials;',targi))
    
    s_z_pool = [];
    s_z_targ_trial_features =  [];
%     for ei = 81
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
%             time_indices = find((s_timePoints <= timeReachStart / 1000) ...
%                 & (s_timePoints >= (timeReachStart - preMovementTime) / 1000)) ;
            freq_indices = find(freq <= freqRange(2)) ;
%             s_z_preMovement = squeeze(mean(s_z(:,freq_indices),2));
            
            %s_z_targ_data(targi,:,:) = s_z_pool_mean;
            s_z_targ_trial_features(ti,ei,:,:) = s_z;
        end
%     end
    
    %build ML datasets
    s_z_targ_ML = cat(1,s_z_targ_ML,s_z_targ_trial_features);
    s_z_targ = [s_z_targ;ones(size(data,1),1) * targi];
    %s_z_targ_features(targi,:,:) = squeeze((s_z_targ_trial_features));
   
    disp(strcat('finished:',num2str(targi)))
    
end
disp('All Processing Done')

%let's remove some redundent dimensions

s_z_targ_ML = squeeze(s_z_targ_ML(:,ei,:,:));

%%
%apply averaging operator
freqStep = 1; % per row
timeStep = 2; % columns
      
s_z_targ_ML_window = [];

for ti = 1:size(s_z_targ_ML,1)
    s_temp = squeeze(s_z_targ_ML(ti,:,:));
    s_temp = s_temp';
    
    s_temp_info.rowIndex = s_timePoints;
    s_temp_info.colIndex = freq(freq_indices);
    
    [B,B_info] = calSlideWindow(s_temp,freqStep,timeStep,s_temp_info);
    s_z_targ_ML_window(ti,:,:) = B';
end

disp('done')

%%

%%
% figure
% histogram(mi_x_y)

figure
imagesc(mi_x_y')
colorbar
axis xy
xlabel('Time from Movement Onset (s)')
ylabel('Frequency (Hz)')
title('Mutual Information between Target and Time Points')

yticklabels = 0:10:300;
yticks = linspace(1, size(mi_x_y, 2), numel(yticklabels));
set(gca, 'yTick', yticks, 'yTickLabel', yticklabels)

xticklabels = B_info.rowIndex(1:8:end) - 2 ;
xticks = linspace(1, size(mi_x_y, 1), numel(xticklabels));
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)

