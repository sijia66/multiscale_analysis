    function trDat = spontMakeTrials(dat, trLen, Fs)
%trDat = spontMakeTrials(dat, trLen, Fs)
%
%divides raw time-series data into arbitrary trials, trLen long
%
%INPUT:  dat   - time-series of data (ch x time)
%        trLen - length of trials (s)
%        Fs    - sampling rate (Hz)
%OUTPUT: trDat - trial-sorted data (ch x time x trials)
%
%A. Orsborn, 2015


%split up into 'trials' of trLen long
nTr = floor(size(dat,2)/ (trLen*Fs));
trDat = nan(size(dat,1), floor(trLen*Fs), nTr); %ch x time x tr

for it = 1:nTr
    ind = (1:floor(trLen*Fs)) + floor(trLen*Fs)*(it-1);
    trDat(:,:,it) = dat(:,ind);
end


