function [trData, eList, t, P_chamber] = loadTrLFP_session(Trials, Session, DriveName, trig, bn, lfpType, monkeydir)

%notes:
%deleted P_drive
global MONKEYDIR

if ~exist('monkeydir', 'var') || isempty(monkeydir)
    monkeydir = MONKEYDIR;
end

if ~exist('lfpType', 'var') || isempty(lfpType)
    lfpType = 'clfp';
end


%get exp def info--pick first recording w/in session. Assumes sessions are grouped by drive-
% --> drive/exp info should the the same across recordings
expDef_file  = [MONKEYDIR '/' Session{1,1} '/' Session{1,2}{1} '/rec' Session{1,2}{1} '.experiment.mat'];
load(expDef_file, 'experiment');
Fs = experiment.hardware.acquisition.samplingrate;

drive_idx = ismember({experiment.hardware.microdrive.name}, DriveName);
eList = [experiment.hardware.microdrive(drive_idx).electrodes.channelid];

%also get exp/electrode position metadata;
%pos = [experiment.hardware.microdrive(drive_idx).electrodes.position];
P_chamber = drmap_coordTransform_microdrive2chamber(experiment, find(drive_idx));

%load trialLfp using info about session
trData = trialLfp(Trials,DriveName,eList,[],trig,bn,lfpType,monkeydir);

%compute time axis for alignment
if Fs==30000
    Fs_lfp = 1000;
else
    Fs_lfp = Fs;
end
t = [bn(1):1e3/Fs_lfp:bn(2)]; 
t = t(1:end-1);
