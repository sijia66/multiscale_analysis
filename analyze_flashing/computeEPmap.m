function [mERP, EPamp, EPampZ, EPmap,EPmapZ,  X, Y, interpMask] = computeEPmap(trData, t, P_chamber, params)

if ~exist('params', 'var') || isempty(params)
    
    params = struct([]);
end

parameters = {'tSearch', 'tBaseline', 'interpMethod', 'interpRadDist', 'eMask', 'rejSigma'}; %rec length in s
defaults   = {[10 100],  [-50 0],      'radialDistance',    2,              'none',   2}; %#ok<NASGU>

fn = fieldnames(params);
for i=1:length(parameters)
    if ~ismember(parameters{i}, fn)
        eval(['params.', parameters{i}, '= defaults{i};'])
    end
end



%%% compute ERPs from trData
mERP = sq(nanmean(trData,1)); %average over trials

%subtract out dc differences across channels
m = sq(nanmean(mERP,2));
mERP = mERP - repmat(m, [1 size(mERP,2)]);


%%% compute ERP amplitude for each electrode
tSearchInds   = t >=params.tSearch(1) & t<=params.tSearch(2);
tBaselineInds = t >=params.tBaseline(1) & t <= params.tBaseline(2);

EP_base = mean(mERP(:,tBaselineInds),2); %baseline = mean over time window for each channel
mxEP    = max(mERP(:,tSearchInds), [],2); %max within search window

EPamp = mxEP - EP_base;



%also z-score EPamps
m = nanmean(EPamp);
sd = nanstd(EPamp);
EPampZ = (EPamp - m)./sd;

%now compute EP map (w/ and w/o z-score)
[Dmap, X, Y] = drmap_getDataMap2(EPamp, P_chamber(1,:), P_chamber(2,:));

[DmapZ, X, Y] = drmap_getDataMap2(EPampZ, P_chamber(1,:), P_chamber(2,:));

%run map interpolation
%set up mask--if specify none, make the mask = ones everywhere
if params.eMask == 'none'
    params.eMask = ones(size(Dmap));
end
%misc specific parameters formatting
if strcmp(params.interpMethod, 'radialDistance')
    interpParams.distRadius = params.interpRadDist;
end
[Dmap_interp, ~] = drmap_spatialMapInterp(Dmap, params.eMask, X, Y, params.interpMethod, interpParams);
[DmapZ_interp, interpMask] = drmap_spatialMapInterp(DmapZ, params.eMask, X, Y, params.interpMethod, interpParams);


%return
EPmap = Dmap_interp;
EPmapZ = DmapZ_interp;
