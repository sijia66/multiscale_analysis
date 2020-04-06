function checkCoherence(day, rec, e1, e2, datType, recLen)


if ~iscell(rec)
    rec = {rec};
end

for r=1:length(rec)
    
    
%load the data up
dat_1 = spontReadTSdat(day, rec{r}, e1.drive, datType, recLen);

if ~strcmp(e1.drive, e2.drive)
    dat_2 = spontReadTSdat(day, rec{r}, e2.drive, datType, recLen);
    
    ind1 = ismember(dat_1.eID,e1.eID);
    ind2 = ismember(dat_2.eID,e2.eID);
    
    driveid1 = ones(1,sum(ind1));
    driveid2 = 2*ones(1,sum(ind2));
    
    dat.tsDat = cat(1, dat_1.tsDat(ind1,:), dat_2.tsDat(ind2,:));
    dat.eID   = cat(2, dat_1.eID(ind1), dat_2.eID(ind2));
    dat.dID   = cat(2, driveid1, driveid2);
    dat.Fs    = dat_1.Fs;
    dat.driveName = cat(2, dat_1.driveName, dat_2.driveName);
    dat.driveNum  = cat(2, dat_1.driveNum, dat_2.driveNum);
    dat.datType   = datType;
else
    
    ind1 = dat_1.eID==e1.eID;
    ind2 = dat_1.eID==e2.eID;
    
    driveid1 = ones(1,sum(ind1));
    driveid2 = 2*ones(1,sum(ind2));
    
    dat.tsDat = cat(2, dat_1.tsDat(ind1,:), dat_1.tsDat(ind2,:));
    dat.eID   = cat(2, dat_1.eID(ind1), dat_1.eID(ind2));
    dat.dID   = cat(2, driveid1, driveid2);
    dat.Fs    = dat_1.Fs;
    dat.driveName = cat(2, dat_1.driveName, dat_1.driveName);
    dat.driveNum  = cat(2, dat_1.driveNum, dat_1.driveNum);
    dat.datType   = datType;
    
end



%generate trials
dat.tsDat_tr = spontMakeTrials(dat.tsDat, 2, dat.Fs);


%use spontCoh to calculate coherence. 
[cohDat(r), ~, ~] = spontCoh(day, rec, [], [], [], dat);
end


%plots for each electrode in drive/group 1, averaged across electrodes in drive/group 2

for i=1:length(e1.eID)
    figure
    hold on
    for r=1:length(rec)
        
        ind1 = dat.dID==1 & dat.eID==e1.eID(i);
        ind2 = dat.dID==2;
        
        plot(cohDat(r).freq, sq( nanmean( real(cohDat(r).coh(ind1,ind2,:)),2)), 'linewidth', 2)
    end
    legend(rec)
    
    set(gca, 'xlim', [0 150], 'ylim', [0 1])
title([e1.drive ' ' num2str(e1.eID(i)) '  ' day])
end



end
    
    



    