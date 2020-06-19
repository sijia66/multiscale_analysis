%% get an example data
% load some example recordings
global MONKEYDIR
%MONKEYDIR = 'D:\SiJia\OneDrive - UW\projects\Brain EEG\data';
MONKEYDIR = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\data'
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

trialInfo.tBefore = -6e2; %ms
trialInfo.tAfter  = 6e2;
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
selectedfile = 'C:\Users\Si Jia\OneDrive - UW\projects\Brain EEG\code\posTargets_180328.mat';
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
%% do bipolar referencing
%use the function multi_rereferencing

P_chamber_ECOG = drmap_coordTransform_microdrive2chamber(experiment,1);
P_chamber_ECOG = P_chamber_ECOG';

P_chamber_SC32 = drmap_coordTransform_microdrive2chamber(experiment,2);
P_chamber_SC32 = P_chamber_SC32';

%only look at good electrodes
P_chamber_ECoG = P_chamber_ECOG(trialInfo.goodECoGs,:);
P_chamber_SC32 = P_chamber_SC32(trialInfo.goodSC32,:);

%plot the relative locations of the electrodes
figure
%first ECoG
scatter(P_chamber_ECoG(:,1), P_chamber_ECoG(:,2),'b')
hold on
scatter(P_chamber_SC32(:,1), P_chamber_SC32(:,2),'r')
hold off
legend('ECoG', 'SC32')
xlabel('mm')
ylabel('mm')
title('Spatial locations of the electrodes')
saveas(gca,'electrode_spatial_locations.png')

%% generate a bipolar referencing scheme
ti = 1; 
distThres = 0.75;

dataForReferencing.data = trLfpData(ti,trialInfo.goodECoGs,:);
dataForReferencing.ePos = P_chamber_ECoG;

%this is a search operation 
resultData = multi_rereferencing...
            (dataForReferencing,distThres);

%plot out te resultant grid
figure
scatter(P_chamber_ECoG(:,1), P_chamber_ECoG(:,2),'b')
hold on
scatter(resultData.ePos(:,1), resultData.ePos(:,2),'r')
legend('ECoG', 'ECoG Multipolar')
xlabel('mm')
ylabel('mm')
title('Spatial locations of the electrodes')
saveas(gca,'electrode_spatial_locations_multipolar.png')


%% 
distThres = 0.75;
[refMatrix_multi, refLocation_multi] = ...
            multi_rereferencing_matrix(P_chamber_ECoG,distThres);
        
P_multi_ref = refLocation_multi * P_chamber_ECoG;

distThres = 0.75;
[refMatrix_bi, refLocation_bi] = ...
            bipolar_matrix(P_chamber_ECoG,distThres);
        
P_bipolar_ref = refLocation_bi * P_chamber_ECoG;

figure
subplot(1,2,1)
scatter(P_chamber_ECoG(:,1), P_chamber_ECoG(:,2),'b')
hold on
scatter(P_multi_ref(:,1), P_multi_ref(:,2),'r')
legend('ECoG', 'ECoG Multipolar')
xlabel('mm')
ylabel('mm')
title('Spatial locations of the electrodes')

subplot(1,2,2)
scatter(P_chamber_ECoG(:,1), P_chamber_ECoG(:,2),'b')
hold on
scatter(P_bipolar_ref(:,1), P_bipolar_ref(:,2),'r')
legend('ECoG', 'ECoG bipolar')
xlabel('mm')
ylabel('mm')
title('Spatial locations of the electrodes')

saveas(gca,'electrode_spatial_multi_bipolar.png')


%%
tca_data = trLfpData(:,trialInfo.goodECoGs,:);
X = tensor(tca_data);
Y = ttm(X,refMatrix_bi,2);
data_temp = double(Y);

%% test the difference between cal_car and matrix
%get only one trial of data
ti = 1;
data_temp_ECoG = squeeze(trLfpData(ti, trialInfo.goodECoGs,:));

%apply car to this single trial of data
%add the car filters
dim = 1;
[data_temp_ECoG_car, data_mean_ECoG]= ...
    cal_car(data_temp_ECoG,dim);


%use the matrix method to get the same data
%we don't squeeze because we are doing a matrix tensor op
data_temp_ECoG_raw = trLfpData(ti, trialInfo.goodECoGs,:);
%get the matrix transformation. 
P_chamber_ECoG = drmap_coordTransform_microdrive2chamber(experiment,2);
P_chamber_ECoG = P_chamber_ECOG(trialInfo.goodECoGs,:);

%get common average referencing filters
refMatrix_car = car_matrix(P_chamber_ECoG);

tca_data = trLfpData(ti,trialInfo.goodECoGs,:);
X = tensor(tca_data);
Y = ttm(X,refMatrix_car,2);
data_temp_tensor = squeeze(double(Y));


data_temp_tensor == data_temp_ECoG_car

figure
imagesc(data_temp_ECoG_car - data_temp_tensor)
colorbar
% subplot(1,2,1)
% imagesc(data_temp_ECoG_car)
% 
% subplot(1,2,2)
% imagesc(data_temp_tensor)


