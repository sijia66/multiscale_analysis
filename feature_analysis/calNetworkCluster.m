function clustercoefficients = calNetworkCluster(cMatrix)
% Input:
%       cMatrix; n by n connectivity matrix
%Output:
%       clustercoefficients connectivty matrix for each electrode
%
% adapted from Analyzing Neural Time series Chapter 31


connmat = cMatrix;

%adopted from analyzing neural time series
for chani=1:size(connmat,1)
    
    % find neighbors (suprathreshold connections)
    neighbors = find(connmat(chani,:));
    n = length(neighbors);
    
    % cluster coefficient not computed for islands
    if n>1
        % "local" network of neighbors
        localnetwork = connmat(neighbors,neighbors);
        % localnetwork is symmetric; remove redundant values by replacing with NaN
        localnetwork = localnetwork + tril(nan(n));
        
        % compute cluster coefficient (neighbor connectivity scaled)
        clustercoefficients(chani,1) = 2*nansum(localnetwork(:)) / ((n-1)*n);
    else
        clustercoefficients(chani,1) = 0;
    end
    
end



end

