function [good_chan] = screenECOGchans(day, rec, params, monkeydir)

global MONKEYDIR

if ~exist('monkeydir', 'var')
    datdir = MONKEYDIR;
elseif isempty(monkeydir)
        datdir = MONKEYDIR;
else
    datdir = monkeydir;
end

if ~exist('rec', 'var')
    REC = dayrecs(day, datdir);
elseif isempty(rec)
    REC = dayrecs(day, datdir);
else
    REC = rec;
end
if ~iscell(REC)
    REC = {REC};
end

if exist('params', 'var')
    fn = fieldnames(params);
else
    fn = {};
    params = struct([]);
end
parameters = {'filetype', 'drive'}; %rec length in s
defaults   = {'lfp', 'all'}; %#ok<NASGU>
for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['params(1).', parameters{i}, '= defaults{i};'])
    end
end



%loop through recordings
for r=1:length(REC)
    
    %load exp def file
    exp_def_name = [datdir '/' day '/' REC{r} '/rec', REC{r}, '.experiment.mat'];
    load(exp_def_name, 'experiment')
    
    drive_names = {experiment.hardware.microdrive(:).name};
    
    if strcmpi( params.drive, 'all' )
        drives = drive_names;
    else
        drives = params.drive;
        
        if ~iscell(drives)
            drives = {drives};
        end
    end
    
    %loop through drives
    for d=1:length(drives)
        
        %load data
        dat = spontReadTSdat(day, REC{r}, drives{d}, params.filetype, [], datdir);
        
        %create trials to compute spectrum
        dat.tsDat_tr = spontMakeTrials(dat.tsDat, 2, dat.Fs);
        
        datParams.fileType = params.filetype;
        datParams.sigmaThreshold = -1; %don't reject outlier trials. too slow for quick inspection.
        %[specDat, ~] = spontPowerSpec([], [], [], 'indChan_trAvg', struct([]), datParams, dat);

        
        %compute correlations
        disp('Calculating correlations...')
        C = spontTScorr(dat.tsDat);
        
        %calculate average correlation for each electrode
        avgC = nanmean(C,2);
        
        %channel ok as long as correlation 2 std fom mean
        up = nanmean(avgC) + 2*nanstd(avgC);
        dn = nanmean(avgC) - 2*nanstd(avgC);
        
        good_chan{d}(r,:) = avgC > dn & avgC < up;
        
        
    end %drives
    
end% recs

  
if length(good_chan)==1
    good_chan = good_chan{:};
end

