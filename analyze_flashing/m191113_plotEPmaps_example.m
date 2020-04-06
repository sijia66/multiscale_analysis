%setting basic info
%clear all


global MONKEYDIR
MONKEYDIR = 'C:/Users/Orsborn Lab/OneDrive - UW/projects/Brain EEG/data';

addpath(genpath([MONKEYDIR '/m/']))

drive_base = 'LM1_ECOG';

%%  get metadata to find days/recs w/ same drive, orientation etc.
%get visEP 'sessions'
taskSessions = VisEP_Database;

%get sessions for 'drive' of interest
driveSessions = makeDriveDatabase(drive_base,{'180328','180328'});

%merge to get sessions w/ labview, split by drive
Sessions = mergeTaskDriveSessions(taskSessions, driveSessions);

%convert to a cell (#sessions x N) rather than nested cells b/c easier  to
%work wiith
Sessions = cat(1, Sessions{:});

%get driveName and theta
driveNames = Sessions(:,3);
driveTheta = [Sessions{:,7}];

%find unique driveNames and the theta for each one
[unique_driveNames, inds] = unique(driveNames);
unique_driveTheta = driveTheta(inds);

nDrives = length(unique_driveNames);


%%

%sub-sampling of sessions if desired
%driveNameAnalyze = unique_driveNames; %all data
driveNameAnalyze = {'LM1_ECOG_3'};
useSess          = ismember(driveNames, driveNameAnalyze);

SessAnalyze       = Sessions(useSess,:);
DayAnalyze        = Sessions(useSess,1);
DriveNamesAnalyze = driveNames(useSess);

%load all relevant trials
nSess = size(SessAnalyze,1);

eTrs = [];

for iD=1:nSess
    %iD
    eTrFN = [MONKEYDIR '/' DayAnalyze{iD} '/mat/ErrorTrials.mat'];
    %size(ErrorTrials)
    load(eTrFN, 'ErrorTrials');
    
    eTrs = cat(2,eTrs, ErrorTrials);
end
ErrorTrials = eTrs;

%visEPs have taskcode = 1 (touch task used)
visStimTrials = ErrorTrials( getTaskCode(ErrorTrials)==1 );

%loop through days to double-check they all have trials
day_list = unique(DayAnalyze);
trialDays = getDay(visStimTrials);

disp(['loaded ' num2str(size(visStimTrials,2)) ' vis stim trials.'])
for iD=1:length(day_list)
    trial_idx = ismember(trialDays, day_list{iD});
    if sum(trial_idx)==0
        day_list{iD}
    end
end
%% now loop through to compute EPs, EP maps, and spline fits for each day 

%trial-alignment parameters
tBefore = -1e3;
tAfter  = 2e3;            
bn = [tBefore tAfter];
trig = 'StartOn';
lfpType = 'lfp';

%EP amp and map parameters
EPparams.tSearch   = [10 100];
EPparams.tBaseline = [-50 0];
EPparams.interpMethod  = 'radialDistance';
EPparams.interpRadDist = 2;
EPparams.eMask         = 'none';
EPparams.rejSigma      = 4;


baseparams.tSearch   = [-300 -150];
baseparams.tBaseline = [-50 0];
baseparams.interpMethod  = 'radialDistance';
baseparams.interpRadDist = 2;
baseparams.eMask         = 'none';
baseparams.rejSigma      = 4;

%parameters for cross-validation folds
%splParams.order = 4;
splParams.order = 6; 
splParams.numKnts = 6;
%trainRatio = 0.4;
%numFolds   = 100;

%loop through days

useTrial = true(length(visStimTrials),1);
for iD=1:length(day_list)
    iD
    disp(['Analyzing day: ' day_list{iD}])
    
    %sub-select trials only belonging to this day
    trial_idx = ismember(trialDays, day_list{iD});
    dayTrials = visStimTrials(trial_idx);
    
    %load all trial data + position info of electrodes
    [trLfpData, eIDs, t, Pchamber_byDay{iD}] = loadTrLFP_session(dayTrials, SessAnalyze(iD,:), DriveNamesAnalyze{iD}, trig, bn, lfpType);
    
    %screen for bad trials
    badTr = screenBadECoGtrials(trLfpData);
    
    nTr_byDay(iD)    = sum(~badTr); %# trials
    nTrBad_byDay(iD) = sum(badTr);
    
    tmp = find(trial_idx);
    useTrial(tmp(badTr)) = false; %also flag trial as bad in master list
     
    %nan-out bad trials across all data
    [nTr, nCh, nT] = size(trLfpData);
    trLfpData(badTr,:,:) = nan;
    
    %also screen for bad electrodes (after removing bad trials)
    badCh = screenBadECoGchannels(trLfpData);
    
    %nan-out bad channels
    trLfpData(:,badCh,:) = nan;
    
    badCh_byDay{iD} = eIDs(badCh);
    nBadCh_byDay(iD) = sum(badCh);
    


%     %normalize trial-lfp data for each electrode into z-scores
%     [a,b,c] = size(trLfpData);
%     tmp = reshape(permute(trLfpData,[2 1 3]), [b a*c]);
%     m = mean(tmp,2)';
%     sd = std(tmp,[],2)';
%     trLfpData = (trLfpData-repmat(m,[a 1 c]))./repmat(sd, [a 1 c]);
%     clear tmp
    
    %compute EP and map in chamber-coordinates
    [mEP_byDay{iD}, EPamp_byDay{iD}, EPampZ_byDay{iD}, EPampMap_chamberCoord_byDay{iD}, EPampMapZ_chamberCoord_byDay{iD}, Xmap_chamber_byDay{iD}, Ymap_chamber_byDay{iD}, interpMask_chamberCoord_byDay{iD}] ...
        = computeEPmap(trLfpData, t, Pchamber_byDay{iD}, EPparams);
    
%     %also compute map in drive-coordinates
%     [~, ~, ~, EPampMap_driveCoord_byDay{iD}, EPampMapZ_driveCoord_byDay{iD}, Xmap_drive_byDay{iD}, Ymap_drive_byDay{iD}, interpMask_driveCoord_byDay{iD}] ...
%         = computeEPmap(trLfpData, t, Pdrive_byDay{iD}, EPparams);
    
    %fit maps with splines
    [splFits_driveCoord(iD,:), EPmapZ_driveCoord_hat(iD,:,:)] = splineApproxMap(Xmap_drive_byDay{iD}, Ymap_drive_byDay{iD}, EPampMapZ_driveCoord_byDay{iD}, splParams);
    
    [splFits_chamberCoord(iD,:), EPmapZ_chamberCoord_hat(iD,:,:)] = splineApproxMap(Xmap_chamber_byDay{iD}, Ymap_chamber_byDay{iD}, EPampMapZ_chamberCoord_byDay{iD}, splParams);
    
    %compute error
    %trim edges b/c of split fit behavior that pins to zero...need to
    %investigate
    [~, normErr_chamberCoord(iD,:), R2_splFit_chamberCoord(iD,:)] = calcMapError(sq(EPmapZ_chamberCoord_hat(iD,2:end-1,2:end-1)), EPampMapZ_chamberCoord_byDay{iD}(2:end-1,2:end-1)); 
    
    [~, normErr_driveCoord(iD,:), R2_splFit_driveCoord(iD,:)] = calcMapError(sq(EPmapZ_driveCoord_hat(iD,2:end-1,2:end-1)), EPampMapZ_driveCoord_byDay{iD}(2:end-1,2:end-1)); 
     
end
%% plot maps

useDay = nTr_byDay>350;
day_list_use = find(useDay);

figure
nPltX = floor(sqrt(length(day_list_use)));
nPltY = nPltX;
if nPltX*nPltY<length(day_list_use)
    nPltX = nPltX+1;
    
    if nPltX*nPltY<length(day_list_use)
        nPltY = nPltY+1;
    end
end
clims_sigma = [-2 2];
cnt=1;
for iD=day_list_use
    
    MAP = EPampMapZ_chamberCoord_byDay{iD}(trimAmount+1:end-trimAmount, trimAmount+1:end-trimAmount) ;
    %MAP = EPampMap_chamberCoord_byDay{iD};
    x = Xmap_chamber_byDay{iD}(trimAmount+1:end-trimAmount);
    y = Ymap_chamber_byDay{iD}(trimAmount+1:end-trimAmount);
    %clims = clims_sigma.*nanstd(MAP(:)) + nanmean(MAP(:));
    
    
    subplot(nPltX, nPltY, cnt)
    %h = imagesc(Xmap_chamber_byDay{iD}, Ymap_chamber_byDay{iD}, MAP, clims);
    h = imagesc(x, y, MAP, [-2 2]);
    %set(h, 'alphadata', 1-interpMask+.8)
    set(gca, 'ydir', 'normal')
    axis square
    title(day_list{iD})
    cnt = cnt+1;
end


%also plot spline fits


figure;

clims_sigma = [-2 2];
cnt=1;
for iD=day_list_use
    
    MAP = sq(splMap_test(iD,:,:));
    %MAP = EPampMap_chamberCoord_byDay{iD};
    
    %clims = clims_sigma.*nanstd(MAP(:)) + nanmean(MAP(:));
    
    
    subplot(nPltX, nPltY, cnt)
    %h = imagesc(Xmap_chamber_byDay{iD}, Ymap_chamber_byDay{iD}, MAP, clims);
    h = imagesc(Xhighres, Yhighres, MAP, [-2 2]);
    %set(h, 'alphadata', 1-interpMask+.8)
    set(gca, 'ydir', 'normal')
    axis square
    title(day_list{iD})
    cnt = cnt+1;
end

