function [specDat, data] = spontPowerSpec(Day, Rec, drive, specType, specParams, datParams, data)

%[specDat, data] = spontPowerSpec(Day, Rec, drive, specType, specParams, datParams);
%
%wrapper function to load and pre-process time-series data and call tfspec to calculate power-spectra 
%loads time-series of data. Divides into arbitrary trials, identifies
%artifact trials, then computes spectra. Options available for how to
%compute spectra. 
%
%INPUTS: Day - string with recording day 'YYMMDD' format
%        Rec - recording number(s) to use. Format options:
%               numeric vector (e.g. [1 2 3])
%               string or cell of strings (e.g. '001' or {'001', '002'})
%               empty or 'all' to use all recordings from a day
%        drive - drives to load. Format options:
%                numeric vector reffering to drive index in exp defn file (e.g. [1 2 3])
%                string or cell of strings with drive names (as set in exp defn file (e.g. 'PMd', or {'PMd', 'M1'})
%        specType - type of spectral estimate to calculate (string). Options:
%                  'indChan_trAvg' - estimate each channel individually, average over trials (created
%                                    arbitrarily by dividing data into chunks) within each channel
%
%        specParams - structure with optional parameters for tfspec.
%           Fields:
%           -'tapers' - tapers to use in [K,TIME], [N,P,K] or [N,W] form.
%			   	        [N,W] Form:  N = duration of analysis window in s.
%                                    W = bandwidth of frequency smoothing in Hz. (default: [1 3 3])
%           -'dn'     - time overlap between neighboring windows
%           -'fk'     - cutoff frequency (Hz) (default: 500)
%           -'pad'    - flag to pad data or not (default: [])
%           -'pval'   - flag to calculate p-values on spectral estimates (default: [])
%        datParams   - structure with optional parameters for loading/pre-processing data        
%           -fileType - type of file to load (e.g. 'raw' , 'lfp', 'clfp'(default)  etc.).
%           -recLen   - length of time to load from each rec file in seconds (default: all)
%           -trLen    - length of 'trials' in seconds (default: 5)
%           -sigmaThreshold - number of SD from mean to use as indicator of an artifact (default: 5)
%                             to SKIP artifact rejection, set sigmaThreshold <= 0
%       data - if data is already loaded, can provide as input to
%              skip re-loading. Expects data to be formatted as outlined in
%              readTSdat. trial-sorted data can be specified with '_tr' tag in the field. 
%OUTPUT:specDat - (#drives x 1) structure. Fields:
%                PS - power spectral data (size depends on specType)
%                freq - frequency values for spectral estimates
%                driveName - name of drive
%                driveNum  - num index of drive (in exp defn file)
%       data - (# drives x 1) structure. Fields:
%               tsDat      - time series (#electrodes x time) data
%               recID      - index of recording file each time-sample belongs to (1 x time)
%               drive      - name of drive
%               Fs         - sampling rate
%               trDat      - time series (#electrodes x time x trials) of trial-sorted data
%               artTr      - flag identifiying artifacts (#electrodes x trials)
%               datType    - fileType of data loaded
%               datParams  - copy of data parameters
%
%A. Orsborn, 2015


if exist('specParams', 'var') 
    if isempty(specParams) || ~isstruct(specParams)
        fn = {};
        specParams = struct([]);
    else
    fn = fieldnames(specParams);
    end
else
    fn = {};
    specParams = struct([]);
end
parameters = {'tapers', 'dn', 'fk', 'pad', 'pval'}; %rec length in s
defaults   = {[1 3 3],  0.05,  500,  [],    []}; %#ok<NASGU>

for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['specParams(1).', parameters{i}, '= defaults{i};'])
    end
end
specParams
if exist('datParams', 'var')
    if isempty(datParams) || ~isstruct(datParams)
        fn = {};
        datParams = struct([]);
    else
        fn = fieldnames(datParams);
    end
else
    fn = {};
    datParams = struct([]);
end
parameters = {'fileType', 'trLen', 'sigmaThreshold', 'recLen'}; %rec length in s
defaults   = {'clfp',      5,        5,                inf}; %#ok<NASGU>

for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['datParams(1).', parameters{i}, '= defaults{i};'])
    end
end

%load up data if not input
if ~exist('data', 'var')
    data = spontReadTSdat(Day, Rec, drive, datParams.fileType, datParams.recLen);
end
nDrive = size(data,2);

specDat(nDrive) = struct('PS', [], 'freq', [], 'driveName', [], 'driveNum', [], 'datType', []);
for iD=1:nDrive
    Fs = data(iD).Fs;
    
    dat = data(iD).tsDat;
    fields = fieldnames(data(iD));
    if all(cellfun('isempty', strfind(fields, '_tr')))
        %divide data into arbitrary trials
        disp('Creating arbitrary trials...')
        trDat = spontMakeTrials(dat, datParams.trLen, Fs);
        data(iD).trSortDat = trDat;
        %eval(['data(iD).', datParams.fileType, '_tr = trDat;'])
    else
        trDat = data(iD).tsDat_tr;
        %eval(['trDat = data(iD).', datParams.fileType, '_tr;'])
    end
    
    if all(cellfun('isempty', strfind(fields, '_artTr'))) & datParams.sigmaThreshold > 0
        %reject artifacts
        disp('Finding outlier trials...')
        artTr = spontFindArtTrials(trDat, datParams.sigmaThreshold);
        data(iD).artTr = artTr;
    elseif ~all(cellfun('isempty', strfind(fields, '_artTr')))
        artTr = data(iD).artTr;
    elseif datParams.sigmaThreshold < 0
        artTr = false(size(trDat,1), size(trDat,3));
    end
    
    data(iD).datParams = datParams;
    disp('computing spectrum...')
    switch lower(specType)
        %compute separately for each channel, averaged across trials (no mean subtraction)
        case {'indchan_travg'}
            
%             hPool = gcp('nocreate');
%             if isempty(hPool)
%                 hPool = parpool('local'); %#ok<NASGU>
%             end
            
              
            for iC = 1:size(trDat,1)
                if mod(iC,10)==1
                    fprintf('     electrode %d...\n', iC)
                end
                %tf spec wants tr/ch x time
                [s, freq] = tfspec(squeeze(trDat(iC,:,~artTr(iC,:)))', specParams.tapers, Fs, specParams.dn, ...
                    specParams.fk,specParams.pad, specParams.pval, 1);
                PS_tmp(iC,:,:) = s;
                freq_tmp(iC,:) = freq;
                
                
            end

            specDat(iD).PS = PS_tmp;
            specDat(iD).freq = freq_tmp(1,:);
    end
    
    specDat(iD).driveName = data(iD).driveName;
    specDat(iD).driveNum  = data(iD).driveNum;
    specDat(iD).datType = datParams.fileType;
end
    
    
    
    
