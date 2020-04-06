function [winX, stdX] = movingWin(x, win, dim, bino)

%[winX, stdX] = movingWin(x, win, dim, bino)
%
%computes a moving average of data x w/ window of win size (in #points)
%If data is binary, estimates probability of true w/ binofit function. 
%Data for 1:win/2 and length(x)-win/2:length(x) is pinned to first
%full-window of data to avoid edge effects. 
%
%inputs: x - data vector. 
%        win - window size (in # points)
%        dim - dimension of x to smooth along (default = 1)
%        bino - flag for binary data (default = false)
%
%outputs: winX - windowed data vector (same size as x)
%         stdX - standard deviation (or 95% CI for binary data) of each window
%
%A. Orsborn 181025

if ~exist('dim', 'var')
    dim = 1;
end
if ~exist('bino', 'var')
    bino = 0;
end

%reshape x to have 'dim' for smoothing as first dimension
nd  = ndims(x);
if nd>6
    warning('movingWin not configured to deal with higher than 6 dimensions. You may want to adjust calculation script!')
end
if nd>2 && bino
    error('movingWin binary option only works for single-dimension vectors')
end

x     = permute(x, [dim, setdiff(1:nd,dim)]);
order = [dim, setdiff(1:nd,dim)];

len = size(x,1);

winX = nan(size(x));
stdX = nan(size(x));
if bino
    stdX = nan(size(x,1),2);
end



for i=ceil(win/2):len-ceil(win/2)
    
    inds = i-ceil(win/2)+1:i+ceil(win/2);
    
    if ~bino
        winX(i,:,:,:,:,:) = squeeze( nanmean(x(inds,:,:,:,:,:),1));
        stdX(i,:,:,:,:,:) = squeeze( nanstd(x(inds,:,:,:,:,:), 1, 1));
    else
        [winX(i), stdX(i,:)] = binofit(sum(x(inds)), win, 0.05);
    end
end

if len>win
    winX(1:ceil(win/2),:,:,:,:,:) = repmat( winX(ceil(win/2),:,:,:,:,:), [ceil(win/2) 1 1 1 1 1]);
    stdX(1:ceil(win/2),:,:,:,:,:) = repmat( stdX(ceil(win/2),:,:,:,:,:), [ceil(win/2) 1 1 1 1 1]);
    
    winX(len-ceil(win/2)+1:end,:,:,:,:,:) = repmat( winX(len-ceil(win/2),:,:,:,:,:), [ceil(win/2) 1 1 1 1]);
    stdX(len-ceil(win/2)+1:end,:,:,:,:,:) = repmat( stdX(len-ceil(win/2),:,:,:,:,:), [ceil(win/2) 1 1 1 1]);
else
    warning('Window length longer than data.')
    size(x)
    win
    winX = repmat( mean(x,1),  [len 1 1 1 1 1]);
    stdX = repmat( std(x,1,1), [len 1 1 1 1 1]);
end

[~, i] = sort(order);

winX = permute(winX, i);
stdX = permute(stdX, i);

    