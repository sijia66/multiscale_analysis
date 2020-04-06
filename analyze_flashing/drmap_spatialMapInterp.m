function [dataMap_interp, interpMask] = drmap_spatialMapInterp(dataMap, eMask, X, Y, interpMethod, interpParams)

% [dataMap_interp, interpMask] = drmap_spatialMapInterp(dataMap, eMask, X, Y, interpMethod, interpParams)
%
%function that takes in a data map (2d matrix) and fills in missing data in
%the map via interpolation. Interpolation is done via averaging nearby
%electrodes, with different options of how to select electrodes to average
%('interpMethod'). It returns both the interpolated data map and a binary
%matrix flagging data that has been interpolated.
%
%input: dataMap - 2d matrix of spatially-organized data. Arranged such that
%                  X-coordinate corresponds to columns, y-coordinate corresponds to rows. 
%       eMask    - binary logical mask (same size as dataMap) to indicate electrodes to be considered in interpolation. 
%                  Use this to exclude electrodes from interpolation calculation. Useful, for instance, to
%                  exclude corerns of 244ch ECOG array that are rounded.
%       X        - vector or 2d matrix w/ X coordinates of each data
%                  element. If a vector, values assumed to correspond to x position
%                  for each column of dataMap.
%       Y        - vector or 2d matrix w/ y coordinates of each data
%                  element. If a vector, values assumed to correspond to y position
%                  for each row of dataMap.
%       interpMethod - string to flag method for selecting nearby electrodes to be averaged for interpolation. 
%                      not case sensitive.
%                      Options:
%                         'nearestneighbor' - closest neighboring electrodes (averages across all w/ same distance).
%                                             If only 1 electrode meets criteria, window is expanded
%                                             to nearest + next-nearest. 
%                         'radialDistance'  - all electrodes w/in radius specified by interpParams.distRadius. 
%                                             distance threshold should be in same units as X, Y
%                                             position variables (typically mm)
%      interpParams - structure for any required parameters for
%                     interpolation methods (see above); optional
%output: dataMap_interp - 2d matrix of spatially-organized data with missing values filled in via interpolation (averaging nearby electrodes)
%        interpMask     - 2d matrix (same size as dataMap) w/ binary
%        logical mask flagging interpolated points (true = interpolated; false = measured).
%
%


%if X and Y are vectors w/ individual lists, replicate into 2d
[ny, nx] = size(dataMap);
if isvector(X)
    %x corresponds to columns in map
    X = X(:)'; %make it a column vector
    
    X = repmat(X, [ny, 1]); %repmat x columns by y rows
end
if isvector(Y)
    %y corresponds to rows in map
    Y = Y(:); %make it a row vector
    
    Y = repmat(Y, [1, nx]); %repmat y rows by x columns
end
    
%find the missing data to be considered (data is a nan + part of mask)
dataMissing = isnan(dataMap) & eMask;
missingData_inds = find(dataMissing);

dataMap_interp = dataMap;

%loop through all missing data points
for iE=missingData_inds(:)'
    
    %electrode position
    epos = [X(iE) Y(iE)];
    
    %cartesian distance to other electrodes
    dist = sqrt( (epos(1)-X).^2 + (epos(2)-Y).^2 );
    dist(iE) = inf; %set current electrode to inf to ignore
    
    
    %a few ways to select electrodes for interpolation
    switch lower(interpMethod)
        
        case 'nearestneighbor'

            %find closest (equidistant)
            nn = find(dist == min(min(dist(~dataMissing & eMask))) & ~dataMissing & eMask);
            
            %if only one, add next-closest
            if length(nn)==1
                dist_tmp = dist;
                dist_tmp(nn) = inf;
                
                nn2 = find(dist_tmp == min(min(dist_tmp(~dataMissing & eMask))) & ~dataMissing & eMask);
                
                nn = cat(1, nn, nn2);
            end
            
            interpEs = nn;
       
        case 'radialdistance'
            
            %find distances < threshold
            interpEs = find( dist < interpParams.distRadius & ~dataMissing & eMask);
            
        otherwise 
            disp('Unsupported interpolation method.')
            return
    end
    
    %replace missing data w/ average of selected electrodes
    dataMap_interp(iE) = mean(dataMap(interpEs));
    
end

%also return a binary mask to flag data that has been interpolated
interpMask = false(size(dataMap));
interpMask(missingData_inds) = true;

            