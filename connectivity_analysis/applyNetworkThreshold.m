function cMatrix = applyNetworkThreshold(cMatrix,threshold)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%let's apply a threshold
cMatrix(isnan(cMatrix)) = 0;
cMatrix(cMatrix > threshold) = 1;
cMatrix(cMatrix ~= 1) = 0;
cMatrix = cMatrix + cMatrix';

end

