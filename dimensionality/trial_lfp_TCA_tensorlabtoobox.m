% load some example recordings
% this does not need to change
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

trialInfo.tBefore = 0; %ms
trialInfo.tAfter  = 100;
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

%% prep TCA analysis

%get rid of the corrupted electrodes 
%either noisy electrodes or lowered electrodes

tca_data = trLfpData;
tca_data(:,trialInfo.badE,:) = [];

%we would also like to re-arrange trials by their reach directions. 



[trial_order_ind, trial_order_val] = ...
                    sort_trial_by_target(trial_TargetAssigned);
                
tca_data = tca_data(trial_order_ind, : ,: );
%we feed this into TCA analysis and try to identify the rank
%we annotate thi data
%num of trials by electrodes by time
% options.MinR = 1;
% options.MaxR = 12;%you have to supply MinR and MaxR at the same time
% rankest(tca_data, options);

%immediate next step would be to zscore all the components


R = 6;
rng(1)

X = tensor(tca_data);
M1 = cp_als(X,R);

vizopts = {'PlotCommands',{'bar','bar','line'},...
    'ModeTitles',{'Trials','Electrodes','Time'},...
    'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
info1 = viz(M1,'Figure',1,vizopts{:});


saveas(gca,'TCA_decomp.png')

%% this is not going anywhere
% can we try only looking at one drive
tca_data = trLfpData(:,trialInfo.goodECoGs,:);
tca_data = tca_data(trial_order_ind, : ,: );

R = 6;
X = tensor(tca_data);
M1 = cp_als(X,R);

vizopts = {'PlotCommands',{'bar','bar','line'},...
    'ModeTitles',{'Trials','Electrodes','Time'},...
    'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
info1 = viz(M1,'Figure',1,vizopts{:});

saveas(gca,['TCA_ECoG_decomp_rank_',num2str(R),'.png'])


%% try this only on SC32

% can we try only looking at one drive
tca_data = trLfpData(:,trialInfo.goodSC32,:);
tca_data = tca_data(trial_order_ind, : ,: );

R = 6;
X = tensor(tca_data);
M1 = cp_als(X,R);

vizopts = {'PlotCommands',{'bar','bar','line'},...
    'ModeTitles',{'Trials','Electrodes','Time'},...
    'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
info1 = viz(M1,'Figure',1,vizopts{:});

saveas(gca,['TCA_SC32_decomp_rank_',num2str(R),'.png'])

%% single electrode analysis

ei_d = 2
%set up parameters
params.Fs = 1000;
tfproduct = 5;
params.tapers = [tfproduct 2*tfproduct-1];
params.trialave = 0;
params.fpass = [0 150];
movingwin = [0.2  0.05];
num_ticks = 4;

dat = squeeze(trLfpData(:,ei_d,:));
data = dat'; % times * trials
%calculate the  spectrum
[S,t,f] = mtspecgramc( data, movingwin, params );
t = t + trialInfo.tBefore/1000; %shift relative to movement start
S_Z = zscore(S,[],1);

%we will first test out the 
tca_data = S_Z;
tca_data = tca_data( : ,:,trial_order_ind );


R = 12;
X = tensor(tca_data);
M1 = cp_als(X,R);

vizopts = {'PlotCommands',{'line','bar','bar'},...
    'ModeTitles',{'Time','frequency','Trials'},...
    'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
info1 = viz(M1,'Figure',1,vizopts{:});

saveas(gca,['TCA_ECoG_tfspe_decomp_electrode_', num2str(ei_d),'_rank_',num2str(R),'.png'])

%we 


