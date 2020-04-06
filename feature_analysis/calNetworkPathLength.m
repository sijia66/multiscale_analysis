function temppathlengths = calNetworkPathLength(connmat)
% Input:
%       cMatrix; n by n connectivity matrix
%Output:
%       pathlength connectivty matrix for each electrode
%
% 
% adapted from Analyzing Neural Time series Chapter 31


temppathlengths = pathlength(double(connmat));
temppathlengths(~isfinite(temppathlengths)) = nan;
temppathlengths = nanmean(temppathlengths,2);

end

