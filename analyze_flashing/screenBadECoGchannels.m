function badCh = screenBadECoGchannels(data)


[a,b,c] = size(data);
if ~any([a b c]==1) %data is 3d (i.e. trial-sorted)
    
    %assume data is time x channels x trials
    %reshape into channels x time
    data = reshape(permute(data,[2 1 3]), [b a*c]);
end
%otherwise, assume data already given in channels x time


%compute standard deviation over time
sd = nanstd(data, [], 2); 
med_sd = nanmedian(sd);

%compute CDF of sd distribution
nbins = size(data,1);
[cnts, bins] = hist(sd, nbins);

CDF = cumsum(cnts)/sum(cnts);

%mark channels with std at top/bottom 5% of distribution as bad
bottomCutoff = bins(find(CDF<=0.05, 1, 'last')); 
topCutoff    = bins(find(CDF>=0.95, 1, 'first'));

badCh = sd<bottomCutoff | sd>topCutoff;

%save channels where the sd is +/- 30% of median
%this avoids unnecessary removal of bad channels if there aren't any
inrange = abs(sd - med_sd)<=5*std(sd(~badCh));

badCh = badCh & ~inrange;

figure
plot(sd, '.')
hold on
plot(find(badCh), sd(badCh), 'r*')
disp(['Flagging ' num2str(sum(badCh)) ' channels as bad'])


