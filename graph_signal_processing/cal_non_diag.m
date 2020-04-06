function vech_off_L= cal_non_diag(vechL,n_nodes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% requires the half vectorization toolbox



M = full(DuplicationM(n_nodes)); % default to lower triangular matrix
L = reshape(M * vechL, n_nodes,n_nodes);
vech_off_L =  L(itril(size(L),-1));



end

