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


%% use this as normalization for an electrode across trials

freqParams.tapers = [0.5 5];
freqParams.Fs = 1000;
freqParams.dn = 0.05;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];


s_timePoints_spacing = 100;
data_all_trials = concat_all_trials(Trials,experiment,driveNameAnalyze,MONKEYDIR);

for ei = 211:211

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
% caxis([-6 6])

xticklabels = s_timePoints(1:s_timePoints_spacing:end);
xticks = 1:s_timePoints_spacing:length(s_timePoints);
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)

%% do self normalization for electrode 211 (TargsOn)
% this run for targsOn
driveNameAnalyze = 'LM1_ECOG_3';
data_all_trials = concat_all_trials(Trials,experiment,'TargsOn',driveNameAnalyze,MONKEYDIR);

freqParams.tapers = [0.5 10];
freqParams.Fs = 1000;
freqParams.dn = 0.05;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];

s_timePoints_spacing = 100;

for ei = 211:211

dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, [], [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);
end

%get the actual data
data_all_trials = concat_all_trials(Trials,experiment,'ReachStop',driveNameAnalyze,MONKEYDIR);
for ei = 211:211

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
title('Normalization Period from StartOn to TargsON for electrode 211')

%% %% do self normalization for electrode 211
% this run for ReachStop change the 

driveNameAnalyze = 'LM1_ECOG_3';
data_all_trials = concat_all_trials(Trials,experiment,'ReachStop',driveNameAnalyze,MONKEYDIR);

freqParams.tapers = [0.5 10];
freqParams.Fs = 1000;
freqParams.dn = 0.05;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];

s_timePoints_spacing = 100;

for ei = 211:211

dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, [], [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);
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
title('Normalization Period from StartOn to ReachStop for electrode 211')


%% write grid_sub 
%% % self normalization and referencing
% our guess is that we are gonna remove all common fluctuations in the
% signal
driveNameAnalyze = 'LM1_ECOG_3';
data_all_trials = concat_all_trials_grid_sub(Trials,experiment,'TargsOn',driveNameAnalyze,MONKEYDIR);

freqParams.tapers = [0.5 10];
freqParams.Fs = 1000;
freqParams.dn = 0.05;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];

s_timePoints_spacing = 100;

for ei = 211:211

dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, [], [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);
end


%get the actual data
data_all_trials = concat_all_trials_grid_sub(Trials,experiment,'ReachStop',driveNameAnalyze,MONKEYDIR);
for ei = 211:211

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
title('Normalization Period from StartOn to TargsON for electrode 211')

%% try the comparision for a different electrode
% a quick note on what what we are doing her.
% we are concatenating from startAcq to TargsOn and use this as
% normalizaton for recording from startAcq to ReachStart

driveNameAnalyze = 'LM1_ECOG_3';
data_all_trials = concat_all_trials_grid_sub(Trials,experiment,'TargsOn',driveNameAnalyze,MONKEYDIR);

freqParams.tapers = [0.5 10];
freqParams.Fs = 1000;
freqParams.dn = 0.1;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];

s_timePoints_spacing = 20;

for ei = 1:211 

dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, freqParams.freqRange, [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);

s_mean_all_trials(ei,:) = s_mean;
s_std_all_trials(ei,:) = s_std;

end
% get the actual data
%get the actual data
data_all_trials = concat_all_trials_grid_sub(Trials,experiment,'ReachStop',driveNameAnalyze,MONKEYDIR);

for ei = 1
    
dat = data_all_trials(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, freqParams.freqRange, [], [], [], freqParams.contFlag);
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
title('Normalization Period from StartOn to TargsON for electrode 1 and Grid Sub')

%looks like we have some periodical band structure in the spectra
%two other factors we have not considered are targets and normalize the
%amplitude in the time domain 


%% try out for 
s_timePoints_spacing = 2;
data = data_ref_mean;
s_z_trials_all_electrodes =  [];

for ei = 1:211
    
    data_electrode = squeeze(data(:,ei,:));
    
    s_mean = squeeze(s_mean_all_trials(ei,:,:));
    s_std = (s_std_all_trials(ei,:,:));
    
    [s_z_trials,s_z_info] = caltf_trialNorm_long(data_electrode,s_mean,s_std,freqParams);
    
    s_z_trials_all_electrodes(ei,:,:,:) = s_z_trials;
  
end

s = squeeze(mean(s_z_trials(:,:,:)));
s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.dn;


figure
imagesc(s')
axis xy
colorbar

xticklabels = s_timePoints(1:s_timePoints_spacing:end);
xticks = 1:s_timePoints_spacing:length(s_timePoints);
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)
title(sprintf('Average Spectrum For Electrode %d',ei))
%%
% we would also need to figure the target and then we can lookat the
% distributions

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

% assignTargetNumbers
trial_TargetAssigned = ...
    assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1);

%% 


%% SC32
% a quick note on what what we are doing her.
% we are concatenating from startAcq to TargsOn and use this as
% normalizaton for recording from startAcq to ReachStart

driveNameAnalyze = 'LM1_SC32_1';
data_all_trials_SC32 = concat_all_trials(Trials,experiment,'TargsOn',driveNameAnalyze,MONKEYDIR);

freqParams.tapers = [0.5 10];
freqParams.Fs = 1000;
freqParams.dn = 0.1;        
freqParams.contFlag = 1;
freqParams.freqRange = [0 300];

s_timePoints_spacing = 20;

for ei = 1:32

dat = data_all_trials_SC32(ei,:);
[s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, ...
                freqParams.dn, freqParams.freqRange, [], [], [], freqParams.contFlag);
s_mean = mean(s);
s_std = std(s);

s_mean_all_trials(ei,:) = s_mean;
s_std_all_trials(ei,:) = s_std;

end
   
%% load the data for SC32 electrodes

s_timePoints_spacing = 2;
data = data_twoDrives{2};
s_z_trials_all_electrodes_SC32 =  [];

for ei = 1:32
    
    data_electrode = squeeze(data(:,ei,:));
    
    s_mean = squeeze(s_mean_all_trials(ei,:,:));
    s_std = (s_std_all_trials(ei,:,:));
    
    [s_z_trials,s_z_info] = caltf_trialNorm_long(data_electrode,s_mean,s_std,freqParams);
    
    s_z_trials_all_electrodes_SC32(ei,:,:,:) = s_z_trials;
    
end

s = squeeze(mean(s_z_trials(:,:,:)));
s_timePoints = freqParams.tapers(1)* 0.5 +(1:s_length)*freqParams.dn;


figure
imagesc(s')
axis xy
colorbar

xticklabels = s_timePoints(1:s_timePoints_spacing:end);
xticks = 1:s_timePoints_spacing:length(s_timePoints);
set(gca, 'xTick', xticks, 'xTickLabel', xticklabels)
title(sprintf('Average Spectrum For Electrode %d',ei))


