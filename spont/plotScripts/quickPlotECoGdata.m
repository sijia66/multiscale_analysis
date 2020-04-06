function quickPlotECoGdata(Day, Rec, drive, datType, recLen, fs)

%revision history
%2019 Jan 24, added a for loop to find the drive in experiment metadata


global MONKEYDIR

if ~exist('recLen', 'var')
    recLen = 30; %length of recording to load in seconds;
end
if ~exist('fMin', 'var')
    fMin = 20;
end
if ~exist('fMax', 'var')
    fMax = 30;
end
if ~exist('fs', 'var')
    fs = 14; %font size for plots
end

disp('Loading data...')
dat = spontReadTSdat(Day, Rec, drive, datType, recLen);

disp('creating trials...')
dat.tsDat_tr = spontMakeTrials(dat.tsDat, 2, dat.Fs);

datParams.fileType = datType;
datParams.sigmaThreshold = -1; %don't reject outlier trials. too slow for quick inspection.
[specDat, ~] = spontPowerSpec([], [], [], 'indChan_trAvg', struct([]), datParams, dat);


figure
subplot(2,2,[1 2])
plot(specDat(1).freq, sq(log(mean(specDat.PS,2)))')
set(gca, 'xlim', [0 200])
xlabel('Frequency (Hz)')
ylabel('Power in Log Scale')
set(gca, 'fontSize', fs)
hold on
plot(specDat(1).freq, sq(log(mean(mean(specDat.PS,2),1))), 'color', [0.1 0.1 0.1], 'lineWidth', 5)



%map electrodes
if isnumeric(Rec)
    tmp = cell(size(Rec));
    for iR=1:length(Rec)
        if Rec(iR)<10
            numTag = '00';
        elseif Rec(iR)<100
            numTag = '0';
        else
            numTag = '';
        end
        
        tmp{iR} = [numTag, num2str(Rec(iR))];
    end
    
    Rec = tmp;
end
load([MONKEYDIR, '/', Day, '/', Rec{:}, '/rec', Rec{:}, '.experiment.mat'], 'experiment');
%experiment.hardware.microdrive(1).pcbConfig = {'3a', '3b', '4a', '4b', '1a', '1b', '2a', '2b'};
% experiment.hardware.microdrive(1).configFile = '/mnt/pesaranlab/People/Amy/ecog/channelMappings/JW_pt244ch';
% [microdrive, row, col, ch_ind, arm, connector] = layoutDef_ECOG244_JWpt(experiment.hardware.microdrive(drive));

%get position data from exp def file
%old drive
%p = [experiment.hardware.microdrive(drive).electrodes(:).position];

%search for the drive number
for di = 1:length(experiment.hardware.microdrive)
    driveTemp = experiment.hardware.microdrive(di).name;
    if(strcmp(driveTemp,drive))
        driveFound = true; %flag found
        driveFoundValue = di;
       break
    end
end
if(~driveFound)
    driveFoundValue = 1;
    f = errordlg('Drive not found in experiment.mat');
end

p = [experiment.hardware.microdrive(driveFoundValue).electrodes(:).position];
row = [p(:).row];
col = [p(:).col];
e_x   = [p(:).x];
e_y   = [p(:).y];

%average over time
P = sq(mean(specDat(1).PS, 2));
f_inds = specDat(1).freq>fMin & specDat(1).freq < fMax;

P_band = sq(mean(log(P(:,f_inds)),2));

%reshape based on row/col
P_band_img = nan(max(row)+1, max(col)+1);
for i=1:length(P_band)
    if ~isnan(row(i))
        P_band_img(row(i)+1, col(i)+1) = P_band(i);
    end
end
subplot(2,2,3)
imagesc(P_band_img)
colorbar
title(['mean log-power in ', num2str(fMin), ' - ', num2str(fMax), ' Hz range'], 'fontSize', fs)
set(gca, 'fontSize', fs)
% 
%corr vs distance
 disp('Calculating correlations...')
 C = spontTScorr(dat.tsDat);
% % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/correlations.csv',C)
% % % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/tsdat.csv',dat.tsDat)
% % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/col.csv',col)
% % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/row.csv',row)
% % % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/microdrive.csv',microdrive)
% % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/ch_ind.csv',ch_ind)
% % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/arm.csv',arm)
% % csvwrite('/mnt/pesaranlab/People/MollerMara/ecog-csv/connector.csv',connector)
dist = computeElectrodeDistance(e_x, e_y);

%reshape into row vector
C = C(:);
dist = dist(:);

nDist = 10;

[~,bins] = hist(dist, nDist);
bins = [0 bins];
for iB=1:nDist
    %all channels
    inds = dist >= bins(iB) & dist < bins(iB+1);
    mC(iB) = nanmean(C(inds),1);
    sdC(iB) = nanstd(C(inds),[],1);
    
end

D = bins(2:end)-diff(bins);

subplot(2,2,4)
plot(D, mC, 'b', 'lineWidth', 2)
hold on
plot(D, mC-sdC, 'b--', 'lineWidth', 2)
plot(D, mC+sdC, 'b--', 'lineWidth', 2)
xlabel('Electrode distance (row/col)', 'fontSize', fs)
ylabel('Mean Correlation', 'fontSize', fs)
set(gca, 'fontSize', fs)

end
