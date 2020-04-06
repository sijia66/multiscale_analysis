function [coh_matrix,f] = cal_coh(data,cohParams)
% this function calculates the coherence matrix for a given set of
% electrodes
%Inputs:
%   data: electrode * time
%   cohParams: inputs to coh
%Outputs
%   coh_matrix: an upper triangular matrix of electrode by electrode
%   f: frequency vector
%note that: this function now returns the complex cross-spectral coherency
%
nE = size(data,1);

for nE1 = 1:nE
    for nE2 = nE1:-1:1
        X  = sq(data(nE1,:));
        Y = sq(data(nE2,:));
        
        [coh,f,~] = ...
        tfcoh(X,Y,cohParams.tapers,cohParams.sampling,cohParams.dn,[],[],[],[],1); % last 1 is for cont signal
        coh_matrix(nE1,nE2,:) = mean(imag(coh),1); %average in time
        
    end
end


end