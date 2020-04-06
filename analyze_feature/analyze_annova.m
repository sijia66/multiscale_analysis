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

trialInfo.tBefore = -1e3; %ms
trialInfo.tAfter  = 1e3;
trialInfo.timeVec = linspace(trialInfo.tBefore,trialInfo.tAfter, trialInfo.tAfter - trialInfo.tBefore);
trialInfo.timeReachStart = abs(trialInfo.tBefore);
trialInfo.trig = 'ReachStart';
trialInfo.lfpType = 'lfp';
trialInfo.Fs_lfp = 1000;
trialInfo.badECoGs = [47 59 163];
trialInfo.ECoG_offset = 211;
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
    trialInfo.goodSC32 = (find( trialInfo.depthProfile > 0))'; 
    trialInfo.goodE = [trialInfo.goodECoGs trialInfo.goodSC32+trialInfo.ECoG_offset];
    trialInfo.badE = setdiff(trialInfo.proc_electrodes,trialInfo.goodE);
    
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

%% build the feature set


disp('test the significant across the electrodes')

example_electrode = 218;
dat = squeeze(trLfpData(:,example_electrode,:));
data = dat'; % times * trials

%set up parameters
params.Fs = 1000;
tfproduct = 5;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 0;
params.fpass = [0 200];
movingwin = [0.2  0.05];

num_ticks = 8;


[S,t,f] = mtspecgramc( data, movingwin, params );
t = t + trialInfo.tBefore/1000; %shift relative to movement start
S_Z = zscore(S,[],1);


figure
imagesc(squeeze(mean(S_Z,3)'))
title('Mean spectrum across all trials')
colorbar
axis xy

add_t_f_axes(ax,t,f,trialInfo,num_ticks)

%% feature selection over all electrodes
disp('Calculate the feature p  values')

feature_p_values_all = [];
%feature_electrode_location = [];
for example_electrode = trialInfo.proc_electrodes 
    
    if ismember(example_electrode,trialInfo.badE)
        %ignore the noisy ECoG and unlowered SC32 electrodes
        feature_p_values_all(example_electrode,:,:) = nan;
        continue
    end

    dat = squeeze(trLfpData(:,example_electrode,:));
    data = dat'; % times * trials
    
    %calculate the  spectrum
    [S,t,f] = mtspecgramc( data, movingwin, params );
    t = t + trialInfo.tBefore/1000; %shift relative to movement start
    S_Z = zscore(S,[],1);
    
    
    feature_p_values = fs_annova(S_Z,trial_TargetAssigned);
    feature_p_values_all(example_electrode,:,:) = feature_p_values;
    
    if mod(example_electrode,10) == 0
        disp(['finished: ',num2str(example_electrode)])
    end
end

save('annova_results_all.mat',...
    'feature_p_values_all')

%% plot out the mean response
 

fpv_ECOG_mean = squeeze(mean(-log10(feature_p_values_all(trialInfo.goodECoGs,:,:))));
fpv_SC32_mean = squeeze(mean(-log10(feature_p_values_all(trialInfo.goodSC32,:,:))));

fpv_ECOG_std = squeeze(std(-log10(feature_p_values_all(trialInfo.goodECoGs,:,:))));
fpv_SC32_std = squeeze(std(-log10(feature_p_values_all(trialInfo.goodSC32,:,:))));


ax = subplot(2,2,1);
imagesc(fpv_ECOG_mean');
colorbar
axis xy
title('mean ECoG p value log first')
add_t_f_axes(ax,t,f,trialInfo,num_ticks)

ax = subplot(2,2,2);
imagesc(fpv_SC32_mean');
colorbar
axis xy
title('mean SC32 p value log first')
add_t_f_axes(ax,t,f,trialInfo,num_ticks)

ax = subplot(2,2,3);
imagesc(fpv_ECOG_std');
colorbar
axis xy
title('std ECoG p value log first')
add_t_f_axes(ax,t,f,trialInfo,num_ticks)

ax = subplot(2,2,4);
imagesc(fpv_SC32_std');
colorbar
axis xy
title('std SC32 p value log first')
add_t_f_axes(ax,t,f,trialInfo,num_ticks)


%% find out the maximum predictive power for a given ROI
feature_p_values_all_log = -log10(feature_p_values_all);

%trialInfo.proc_electrodes
t_ROI = [-0.5 0.5];

t_ind = (t >= t_ROI(1)) & (t <= t_ROI(2));
t_ROI_arr = t(t_ind);

%initialize the arrays
feature_electrode_max = [];
feature_electrode_location = [];

for example_electrode = trialInfo.proc_electrodes
    f_p_log = squeeze(feature_p_values_all_log(example_electrode,t_ind,:));
    
    if ismember(example_electrode,trialInfo.badE)
        %ignore the noisy ECoG and unlowered SC32 electrodes
        feature_electrode_max(example_electrode) = nan;
        
        feature_electrode_location(example_electrode,:)...
            = [nan nan];
        
        continue
    end
    
    [M,index] = max(f_p_log,[],'all','linear');
    [row,col] = ind2sub(size(f_p_log),index);
    
    feature_electrode_max(example_electrode) = M;
    
    feature_electrode_location(example_electrode,:)...
                            = [t_ROI_arr(row) f(col)];
    
    %[M_list,Index] = sort(f_p_log(:),'descend');
    
end

%%
% examine the data
figure
scatter(trialInfo.proc_electrodes, feature_electrode_max)
title('Max log transformed p-values')

figure
subplot(2,1,1)
scatter(trialInfo.proc_electrodes, feature_electrode_location(:,1))
title('Time at which max predictive value occurs')
ylabel('Time (s)')

subplot(2,1,2)
scatter(trialInfo.proc_electrodes, feature_electrode_location(:,2))
title('Frequency at which max predictive value occurs')
ylabel('freqency (Hz)')
xlabel('Electrode')

saveas(gca, 'feature_map_time_frequency.png')


%% report results in a matlab figure
%get electrode info

P_chamber_ECOG = drmap_coordTransform_microdrive2chamber(experiment,1);
P_chamber_ECOG = P_chamber_ECOG';

P_chamber_SC32 = drmap_coordTransform_microdrive2chamber(experiment,2);
P_chamber_SC32 = P_chamber_SC32';

%looking at significance of the time points
[map_ECoG,mapInfo] = scatterToMatrix(P_chamber_ECOG,...
                        feature_electrode_location(trialInfo.ECOG_indices,1));

figure('units','normalized','outerposition',[0 0 1 1])
set(gca,'color',0*[1 1 1]);

fig = subplot(2,2,1);        
imagesc(map_ECoG', 'AlphaData',double(~isnan(map_ECoG')))
title('Time at which target separation is at maximum for ECoG ')
colormap(redblue)
colorbar
set(gca,'color',0*[1 1 1]);
scatterToMatrix_plot_label(fig, map_ECoG, mapInfo)

%looking at significance of frequency bands
[map_ECoG,mapInfo] = scatterToMatrix(P_chamber_ECOG,...
                        feature_electrode_location(trialInfo.ECOG_indices,2));
fig = subplot(2,2,2);        
imagesc(map_ECoG', 'AlphaData',double(~isnan(map_ECoG')))
title('Freq at which target separation is at maximum for ECoG ')
colormap(redblue)
colorbar
set(gca,'color',0*[1 1 1]);
scatterToMatrix_plot_label(fig, map_ECoG, mapInfo)


%looking at SC32 electrode
% frequency 

[map_SC32,mapInfo] = scatterToMatrix(P_chamber_SC32,...
                        feature_electrode_location(trialInfo.SC32_indices,1));
fig = subplot(2,2,3);       
imagesc(map_SC32', 'AlphaData',double(~isnan(map_SC32')))
title('Time at which target separation is at maximum for SC32 ')
colormap(redblue)
colorbar
set(gca,'color',0*[1 1 1]);
scatterToMatrix_plot_label(fig, map_SC32, mapInfo)

%freq
fig = subplot(2,2,4);    
[map_SC32,mapInfo] = scatterToMatrix(P_chamber_SC32,...
                        feature_electrode_location(trialInfo.SC32_indices,2));
imagesc(map_SC32', 'AlphaData',double(~isnan(map_SC32')))
title('Freq at which target separation is at maximum for SC32 ')
colormap(redblue)
colorbar
set(gca,'color',0*[1 1 1]);
scatterToMatrix_plot_label(fig, map_SC32, mapInfo)

saveas(gca, 'feature_spatial_time_frequency.png')


%% add behaviour info to this data

%plot_trials(Trials)
%make a list of xlines to plot 

%units in ms
startAcq = [Trials.StartAq] - [Trials.ReachStart];
targsOn = [Trials.TargsOn] - [Trials.ReachStart];
go = [Trials.Go] - [Trials.ReachStart];
reachStop = [Trials.ReachStop] - [Trials.ReachStart];

%insert zero for reachStart
behav_times = [startAcq' go' targsOn' zeros(length(Trials),1) reachStop'];

behav_times_mean = nanmean(behav_times) / 1000; %convert to in s
behav_times_std = nanstd(behav_times) / 1000;   %convert to in s

behav_struct.times = behav_times_mean;
behav_struct.times_std = behav_times_std;
behav_struct.labels = {'Start acq' 'Go' 'Targ On' 'Reach Start' 'Reach Stop'};

%% example electrodes at different times
example_electrodes = [1,11,218,211+13];


%set up parameters
params.Fs = 1000;
tfproduct = 5;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 1;
params.fpass = [0 200]; %controls the output frequency window
movingwin = [0.2  0.05];
num_ticks = 4;




%make a new figure
figure('units','normalized','outerposition',[0 0 0.5 1])


for ei = 1:length(example_electrodes)
    ei_d = example_electrodes(ei);

    dat = squeeze(trLfpData(:,ei_d,:));
    data = dat'; % times * trials
    %calculate the  spectrum
    [S,t,f] = mtspecgramc( data, movingwin, params );
    t = t + trialInfo.tBefore/1000; %shift relative to movement start
    S_Z = zscore(S,[],1);
    
    feature_p_values = squeeze(-log10(feature_p_values_all(ei_d,:,:)));
    
    t_window = [behav_struct.times(1) behav_struct.times(end)];
    t_index  = find((t >= t_window(1)) & (t <= t_window(2)));
    
    subplot(4,2,2*ei - 1)
    imagesc(S_Z(t_index,:,:)')
    if ei_d > 211
        ei_d = ei_d - 211;
        title(['Mean spectrum across all trials for SC32 ', num2str(ei_d)])
    else
        title(['Mean spectrum across all trials for ECoG ', num2str(ei_d)])
    end
    colorbar
    axis xy
    %add_t_f_axes(ax,t,f,trialInfo,num_ticks)
    add_t_f_axes_behavInfo(fig,t(t_index),f,behav_times_mean,trialInfo,num_ticks)
    
    %plot out the feature p_values
    %3 is when is target on
    t_window = [behav_struct.times(3) behav_struct.times(end)];
    t_index  = find((t >= t_window(1)) & (t <= t_window(2)));
    
    subplot(4,2,2*ei)
    imagesc( feature_p_values(t_index,:)')
    if ei_d > 211
        ei_d = ei_d - 211;
        title(['Log p values  for SC32 ', num2str(ei_d)])
    else
        title(['Log p values  for ECoG ', num2str(ei_d)])
    end
    colorbar
    axis xy
    %add_t_f_axes(ax,t,f,trialInfo,num_ticks)
    add_t_f_axes_behavInfo(fig,t(t_index),f,behav_times_mean(3:end),trialInfo,num_ticks)
        
end
    
saveas(gca, 'feature_example_time_frequency.png')





