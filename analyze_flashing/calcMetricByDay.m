function [M, SD, day_list] = calcMetricByDay(dat, dayid, metric, n, seg)


if ~exist('metric', 'var') || isempty(metric)
    metric = 'mean';
end

if ~exist('n', 'var') || isempty(n)
    n = 'all';
end
if ~exist('seg', 'var') || isempty(seg)
    seg = 'first';
else
    if ~ismember(seg, {'first', 'last'})
        error('seg can only be first or last')
    end
end


%
day_list = unique(dayid);
nD = length(unique(dayid));

M = nan(nD,1);
SD = nan(nD,1);

%loop over days
for d=1:nD
    
    %indices belonging to given day
    inds = dayid==day_list(d);
    
    %get portion of data to calc over
    if ~strcmp(n, 'all')
        inds = find(inds==1, n, seg);
    else
        inds = find(inds==1);
    end
    
    switch metric
        case 'mean'
            M(d) = nanmean(dat(inds));
            SD(d) = nanstd(dat(inds))./sqrt(numel(inds));
            
        case 'median'
            M(d) = median(dat(inds));
            SD(d) = std(dat(inds))./sqrt(numel(inds));
        case 'max'
            M(d) = max(dat(inds));
            SD(d) = [];
        case 'min'
            M(d) = min(dat(inds));
            SD(d) = [];
        case 'bino'
            try
                [M(d), SD(d)] = binofit(sum(dat(inds)), numel(inds), 0.05);
            catch
                M(d) = binofit(sum(dat(inds)), numel(inds), 0.05);
            end
        otherwise
            disp('Selected metric not supported by calcMetricByDay.')
    end
    
    
end
  


    
    
    