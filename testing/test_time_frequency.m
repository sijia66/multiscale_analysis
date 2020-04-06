%%
% load some example recordings
global MONKEYDIR
MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
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

trialInfo.tBefore = -2e3;
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;
trialInfo.badECoGs = [47 59 163];
trialInfo.goodECoGs = setdiff(1:211,trialInfo.badECoGs);


%we are going to look at all electrode,  211 ECOG and 32 SC32 electrodes
trialInfo.proc_electrodes = 1:243;

%use constant target defination
selectedfile = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
load(selectedfile)

%we load the regular trials
SessAnalyze = Sessions(useSess,:);
nSess = size(SessAnalyze,1);

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
end


%% first, we test the multitaper methods
% some parameters to vary
% Amy's PhD thesis:  200 ms sliding window with a 50ms step-size,
% time-bandwidth product of 3 and 5 tapers, z scored
% 
% Watanabe et al reconstruction: non-overlapping sliding window of 300 ms

example_electrode = 1;
dat =  squeeze(trLfpData(:,example_electrode,:));
disp(['Using all data for electrode:',num2str(example_electrode)])

tapers_time_range = 0.1:0.1:0.5;
tapers_freq_range = 10:10:50;

N_TIMES = length(tapers_time_range);
N_FREQ = length(tapers_freq_range);

%plotting setup
num_ticks = 4


%Amy's method

freqParams.Fs = 1000;
freqParams.freq = 200; %Hz
freqParam.dn = 0.05;
freqParams.pad = 4;
freqParams.flag = 0; %if pooling across all trials 
freqParams.contFlag = 1;

figure('units','normalized','outerposition',[0 0 1 1])
plot_counter = 0;
for i = 1:N_TIMES
    for j = 1:N_FREQ
        
        %parameter sweep
        freqParams.tapers = [tapers_time_range(i) tapers_freq_range(j)];
        
        %calculate the spectrum
        [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, freqParams.freq, freqParams.pad,    [],  freqParams.flag    ,freqParams.contFlag,[]);
        %tfspec(X,             tapers,      sampling,           dn,              fk,             pad, pval,   flag           ,            contflag, Errorbar)
        % do z_scoring
        s_t = reshape(s,size(s,2) * size(s,1),size(s,3));
        s_z_score = zscore(s_t,[],1);
        s_z_mean = sq(mean(reshape(s_z_score,size(s,1), size(s,2) ,size(s,3))));
        
        %s_z_all(i,j,:,:) = s_z_mean;
        
        %then plot
        plot_counter =  plot_counter + 1;
        subplot(N_TIMES,N_FREQ,plot_counter)
        imagesc(s_z_mean')
        axis xy
        title([num2str(tapers_time_range(i)),' s ',... 
               num2str(tapers_freq_range(j)), ' Hz z',...
               't-f prod ', num2str(tapers_time_range(i) * tapers_freq_range(j))])
        

        
        %add time points 
        [t_length, f_length] = size(s_z_mean);
        %add movement start
        xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);
        
        s_timePoints = freqParams.tapers(1)* 0.5 +(1:t_length)*freqParam.dn;
        s_timePoints = s_timePoints + trialInfo.tBefore / 1000; 
        
        tick_step = round(t_length / num_ticks);
        xticks(1:tick_step:t_length)
        xticklabels(s_timePoints(1:tick_step:end))
        
        tick_step = round(f_length / num_ticks);
        yticks(1:tick_step:f_length)
        yticklabels(round(freq(1:tick_step:f_length)))
        
        
    end
    disp(['Finished: ',num2str(i/length(tapers_time_range))])
end

dim = [0 0 .1 .1];
str = ['electrode:',num2str(example_electrode),' mean zscored all trials'];
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'FontSize', 9);

saveas(gca,'pesaran_time_frequency_sweep.png')

%% overlapping sweep step size
disp('overlapping step size')
%fill out the best parameters for the tapers 
freqParams.tapers = [0.2 25];

tn_range = [0.02 0.04 0.05 0.1 0.2 0.5]; 
N_TN = length(tn_range);

freqParams.Fs = 1000;
freqParams.freq = [0 200]; %Hz
freqParams.pad = 4;
freqParams.flag = 0; %if pooling across all trials 
freqParams.contFlag = 1;

figure('units','normalized','outerposition',[0 0 1 0.25])
for i = 1:N_TN
    freqParam.dn = tn_range(i);
    
    %calculate the spectrum
    [s, freq] = tfspec(dat, freqParams.tapers, freqParams.Fs, freqParam.dn, freqParams.freq, freqParams.pad,    [],  freqParams.flag    ,freqParams.contFlag,[]);
    %tfspec(X,             tapers,      sampling,           dn,              fk,             pad, pval,   flag           ,            contflag, Errorbar)
    
    % do z_scoring
%     s_t = reshape(s,size(s,2) * size(s,1),size(s,3));
%     s_z_score = zscore(s_t,[],1);
%     s_z_mean = sq(mean(reshape(s_z_score,size(s,1), size(s,2) ,size(s,3))));

    %do z-scoring
    s_z = zscore(s,[],2);
    s_z_mean = squeeze(mean(s_z));
    
    
    %then plot
    subplot(1,N_TN,i)
    imagesc(s_z_mean')
    axis xy
    title({
        [num2str(freqParams.tapers(1)),' s ',...
        num2str(freqParams.tapers(2)), ' Hz z',...
        't-f prod ', num2str(prod(freqParams.tapers))];...
        ['step size ',num2str(freqParam.dn),' s']...
        })
    
    %add time points
    [t_length, f_length] = size(s_z_mean);
    %add movement start
    xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);
    
    s_timePoints = freqParams.tapers(1)* 0.5 +(1:t_length)*freqParam.dn;
    s_timePoints = s_timePoints + trialInfo.tBefore / 1000;
    
    tick_step = round(t_length / num_ticks);
    xticks(1:tick_step:t_length)
    xticklabels(s_timePoints(1:tick_step:end))
    
    tick_step = round(f_length / num_ticks);
    yticks(1:tick_step:f_length)
    yticklabels(round(freq(1:tick_step:f_length)))
    
end


saveas(gca,'time_step_sweep.png')

%% test the methods from chronux
disp('Remember to set the path for this toolbox')

taper_range = [1 3 5 10 20 25];
NUM_TAPERS = length(taper_range);
disp('Taper in the form of [NW 2NW-1]')

dat = squeeze(trLfpData(:,example_electrode,:));
data = dat'; % times * trials

params.Fs = 1000;
params.trialave = 1;
params.fpass = [0 200];

movingwin = [0.2  0.05];

figure('units','normalized','outerposition',[0 0 1 0.25])

for i = 1:NUM_TAPERS
    
    params.tapers = [taper_range(i) (taper_range(i) * 2-1)];
    
    [S,t,f] = mtspecgramc( data, movingwin, params );
    %S is the form of time * freq * trials
    
    S = zscore(S,[],1);
    subplot(1,NUM_TAPERS,i)
    imagesc((S'))
    axis xy
    title(['NW: ', num2str(taper_range(i))])
    
    %add time points
    t_length = length(t);
    f_length = length(f);
    
    %add movement start
    xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);
  
    s_timePoints = t;
    freq = f;
    s_timePoints = s_timePoints + trialInfo.tBefore / 1000;
    
    tick_step = round(t_length / num_ticks);
    xticks(1:tick_step:t_length)
    xticklabels(s_timePoints(1:tick_step:end))
    
    tick_step = round(f_length / num_ticks);
    yticks(1:tick_step:f_length)
    yticklabels(round(freq(1:tick_step:f_length)))
end
saveas(gca,'chronux_tfprod_sweep.png')
disp('Remember to delete the paths for this toolbox')

%% within multitaper, we can test overlapping or non-overlapping windows

%% another way to do this is to do frequency filtering