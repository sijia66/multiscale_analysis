function badTr = screenBadECoGtrials(data, sdThresh, chFrac)

% badTr = screenBadECoGtrials(data, sdThresh, chFrac)
%screens for bad trials in trial-aligned continuous data
%computes the mean and std of each channel across all time (tr x time), and
%finds time-points that are sdThresh*std away from the mean. A trial is
%flagged as bad for a given channel if any time-pom

if ~exist('sdThresh', 'var') | isempty(sdThresh)
    sdThresh = 5;
end
if ~exist('chThresh', 'var') | isempty(chFrac)
    chFrac = .5;
end
%data = [trials x channels x time]
[nTr,nCh,nT] = size(data);

data = reshape(permute(data, [2 1 3]), [nCh nTr*nT]);


%compute median and standard deviation over time
sd  = nanstd(data, [], 2); 
med = nanmedian(data, 2);

%find time-points outside of range
badTimePoints = abs(data - repmat(med, [1 nTr*nT]))>repmat(sdThresh.*sd, [1 nTr*nT]);

%reshape into trials x channels x time
badTimePoints = permute( reshape(badTimePoints, [nCh nTr nT]), [2 1 3]);

%find trials w/ any bad time-points across all channels
badTrialbyCh = any(badTimePoints, 3);

badTr = sum(badTrialbyCh,2)>chFrac*nCh;

disp(['Flagging ' num2str(sum(badTr)) ' trials of ' num2str(sum(nTr)) ' as bad'])


