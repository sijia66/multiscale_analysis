function [filtered_data, data_mean]= cal_car(data,dim)
% this function subtracts out the common average across all electrodes
% Input:
%       data: num of trials by num_electrodes by num_time
%       dim; which dimension to subtract out
% Output:
%       filterd_data: num_trials by num_electrodes by num_time
%       data_mean: num_trials by num_electrodes by num_time
% example usage:     
%       trLfpData_ECOG_CAR = cal_car(trLfpData(:,trialInfo.ECOG_indices,:),2)

data_size = size(data);
data_size_logic = (data_size >= 1).*1;
data_size_logic(dim) = data_size(dim);

data_mean = mean(data, dim);
filtered_data = data -  repmat(data_mean,data_size_logic);

end