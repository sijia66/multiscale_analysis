function [mEP, mEP_artrej, dat, dat_artrej, t, pulseLabels] = computeEvokedPotentials(day, rec, drive, params, M, SD)

global MONKEYDIR
Fs_raw = 3e4; %sampling rate for raw files, comedi files

if exist('params', 'var')
    if ~isempty(params)
        fn = fieldnames(params);
    else
        fn = {};
    end
else
    fn = {};
    params = struct([]);
end
parameters = {'bn',      'datType', 'triggerChannel', 'sigmaThreshold', 'normalize'}; %rec length in s
defaults   = {[-20 100],  'lfp',     4,                4,                 1}; %#ok<NASGU>

for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['params(1).', parameters{i}, '= defaults{i};'])
    end
end



datdir = [MONKEYDIR '/' day '/' rec '/'];
expDefname = [datdir 'rec', rec, '.experiment.mat'];

load(expDefname, 'experiment');
if isnumeric(drive)
    drivename = experiment.hardware.microdrive(drive).name;
    drivenum = drive;
elseif ischar(drive)
    drivename = drive;
    tmp = {experiment.hardware.microdrive(:).name};
    drivenum =  find(ismember(tmp, drive));
else
    error('drive must be a name or number')
end

nCh = length(experiment.hardware.microdrive(drivenum).electrodes);        %num channels
%ePos = [experiment.hardware.microdrive(drivenum).electrodes(:).position]; %position data


datfile = [datdir, 'rec', rec, '.', drivename, '.' params.datType, '.dat'];
comedifile = [datdir, 'rec', rec, '.comedi.dat'];
switch params.datType
    case 'raw'
        Fs = 3e4;
        datFormat = 'short';
    case {'lfp', 'clfp', 'mlfp'}
        Fs = 1e3;
        datFormat = 'float';
end

dat_fid = fopen(datfile);
data    = fread(dat_fid, [nCh inf], datFormat);
fclose(dat_fid);

%if normalizing data, do so
if params.normalize
    data = rms_norm_dat(data);
end

com_fid = fopen(comedifile);
comedi_dat = fread(com_fid, [8 inf], 'ushort');
fclose(com_fid);

%get trigger data, adjust sampling rate etc.
trig_dat = comedi_dat(params.triggerChannel,:);
trig_dat = trig_dat(1:Fs_raw/Fs:end);
%trig_dat = smooth(trig_dat,10);


%get evoked potentials
[st, rmvTr] = getEvokedPotentials(data, trig_dat, 1, params.bn, Fs, []);

%normalize data (zscore) -- if not provided, unit transform
if ~exist('M', 'var') || isempty(M)
    M = zeros(nCh,1);
else
    M = M(:);
end

if ~exist('SD', 'var') || isempty(SD)
    SD = ones(nCh,1);
else
    SD = SD(:);
end
st = (st - repmat(M, [1 size(st,2) size(st,3)]))./ repmat(SD, [1 size(st,2) size(st,3)]);

%also label pulse times
pulseLabels = labelPulseTypes(trig_dat, Fs);
pulseLabels = pulseLabels(~rmvTr);

%mean EP
mEP = mean(st,3);

%zero potentials @ start of trial before taking average
t = params.bn(1):1000/Fs:params.bn(2);
ind = find(t<0, 1, 'last');

st = st - repmat( st(:,ind,:), [1 size(st,2), 1]);
mEP_norm = mean( st, 3);

%reject 'artifacts'
st_corr = removeArtifacts_trialBased(st, params.sigmaThreshold);

%mean
mEP_stim_artrej = nanmean(st_corr,3);

%zero @ start
x = st_corr - repmat( st_corr(:,ind,:), [1 size(st_corr,2), 1]);
mEP_norm_artrej = nanmean( x, 3);



mEP = mEP_norm; %_artrej;
dat = st; %_corr;

mEP_artrej = mEP_norm_artrej;
dat_artrej = st_corr;

