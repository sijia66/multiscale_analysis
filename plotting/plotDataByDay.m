function [axH, figH, M, SD] = plotDataByDay(dat, dayid, sm, bino, axH, colors, ylim, plotFlag)

% [axH, figH] = plotMetricByDay(dat, dayid, sm, bino, axH, colors)
%
%This function takes a data vector and IDs for day the data belongs to and
%generates a plot of the (smoothed) data separated by day.
%You can use this, for instance, to plot the evolution of reach time across
%trials within and across behavioral sessions. 
%
%example use flow: 
%1)load & concatonate Trials structure across all day(s).
%2) Calculate reach time w/ a calc function operating on the Trials
%3) find dayid for each Trial in structure
%4) plotDataByDay to make a plot of reach time over days. 
%
%inputs: dat   - (# trials x 1) data vector to plot
%        dayid - (# trials x 1) vector with ID of day each data point belongs to
%        sm    - smoothing window length (# trials). Set to empty to plot
%                raw, unsmoothed data. Smoothing is done WITHIN day to avoid
%                blurring effects across days. 
%        bino  - flag for binary-variable data (e.g. whether a trial is
%                successful or not). Binary data is plotted as probability of true,
%                estimated with matlab 'binofit' function. Default: 0
%        axH   - handle to an existing figure axis. Use this input to have
%                function generate the plot on an existing figure. If not
%                specified, function generates a new figure and returns the handle.
%        colors - (#colors x 3) matrix of RGB color values to use for the
%                 plot (one color / day). #colors must be >= #days in data. 
%                 default: cbrewer Green -> blue color progression, len =
%                 10
%        plotFlag - flag to turn off plotting if desired (default: 1 (make plot))
%outputs: axH  - handle of axis for plot. 
%         figH - handle of figure for plot.
%         M    - mean data (smoothed data plotted)
%         SD   - standard dev (or CI) of plotted data
%
%
%A. Orsborn, 181025

if ~exist('plotFlag', 'var')
    plotFlag = 1;
end

if exist('axH', 'var') && ~isempty(axH) && plotFlag
    figH = get(axH, 'Parent');
    figure(figH)
    set(figH,'CurrentAxes', axH);
    hold(axH)
    
    figH = get(axH, 'parent');
elseif plotFlag
    figH = figure;
    axis();
    axH = get(figH, 'Children');
    hold(axH)
end
if ~exist('colors', 'var') || isempty(colors)
    colors = cbrewer('seq', 'GnBu', 10);
    colors = colors(3:end,:);
end
if exist('sm', 'var') && ~isempty(sm)
    smoothFlag = 1;
else
    smoothFlag = 0;
end

if ~exist('bino', 'var')
    bino = 0;
end

%
day_list = unique(dayid);
nD = length(unique(dayid));

%get range of data for plotting purposes
if smoothFlag
    %tmp = smooth(dat,sm);
    tmp = movingWin(dat, sm,1, bino);
else
    tmp = dat;
end
if exist('ylim', 'var')
    plt_mn = ylim(1);
    plt_mx = ylim(2);
else
    mx = max(tmp); mn = min(tmp);
    if mn<0
        plt_mn = mn*2;
    else
        plt_mn = mn*.5;
    end
    if mx>0
        plt_mx = mx*2;
    else
        plt_mx = mx*.5;
    end
end

M = nan(length(dat),1);
if bino
    SD = nan(length(dat),2);
else
    SD = nan(length(dat),1);
end

%loop over days
for d=1:nD
    
    %indices belonging to given day
    inds = find(dayid==day_list(d));
    
    %grab 
    if smoothFlag
        [x, sdx] = movingWin(dat(inds), sm, 1, bino);
        %x = smooth(dat(inds), sm);
        
        M(inds,:) = x;
        SD(inds,:) = sdx;
    else
        x = dat(inds);
        sdx = [];
        
        M(inds,:) = x;
    end
    
    
    
    if plotFlag
        plot(axH, inds, x, 'color', colors(d,:), 'linewidth', 3)
        if ~isempty(sdx)
            if size(sdx,2)==1
                ciplot(x-sdx./sqrt(sm), x+sdx./sqrt(sm), inds, colors(day_list(d),:))
            elseif size(sdx,2)==2
                ciplot(sdx(:,1), sdx(:,2), inds, colors(day_list(d),:))
            end
        end
        
        %divider among days
        if d<nD
            plot(axH, [inds(end) inds(end)], [plt_mn plt_mx], 'k--')
        end
    end
end
   
if plotFlag
    %set axes tight around max/min OR around user-specified ylim
    if exist('ylim', 'var')
        plt_mn = ylim(1);
        plt_mx = ylim(2);
    else
        
        if mn<0
            plt_mn = mn*1.1;
        else
            plt_mn = mn*.9;
        end
        if mx>0
            plt_mx = mx*1.1;
        else
            plt_mx = mx*.9;
        end
    end
    set(axH, 'ylim', [plt_mn plt_mx], 'xlim', [1 length(dayid)])
end


    
    
    