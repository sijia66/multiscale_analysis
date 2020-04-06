%combine array
function data_out = compress_dim(data,dim)
% this function squashes the defined dimension and then moves the trial col
% by colum
%input:
%   data: 3D matrix
%   dim: dimension to be preserved
%
%example
%data = SC32_lowered_lfp;

z = size(data);

zl = 1:length(z);
zl(dim) = [];

data_temp = permute(data, [zl dim]);


for ei = 1:z(dim)
    data_i = squeeze(data_temp(:,:,ei));
    data_i = reshape(data_i',[],1);
    
    data_out(:,ei) = data_i;
    
end

end