function data_ref_mean = ref_grid_sub(data_temp)
% this function subtracts grid average from a single trial recording.
% note that the grid is hardcoded
% Inputs:
%       data_temp: recordings in the form [electrode time]
% Output:
%       data_ref_mean: grid subtracted data in [electrode time]

% ti = 1:length(Trials); 
outlierElectrodes = [47,59,163];
% fprintf('Set the following electrodes to NAN:'); outlierElectrodes

data_temp(outlierElectrodes,:) = nan;

gridCoor_1 = 1:61;
gridCoor_2 = 62:90;
gridCoor_3 = 91:151;
gridCoor_4 = 152:211;


dataForReferencing(1).data = (data_temp(gridCoor_1,:));
dataForReferencing(1).mean = nanmean(dataForReferencing(1).data,1);
dataForReferencing(1).data_ref_mean = dataForReferencing(1).data - dataForReferencing(1).mean;

dataForReferencing(2).data = (data_temp(gridCoor_2,:));
dataForReferencing(2).mean = nanmean(dataForReferencing(2).data,1);
dataForReferencing(2).data_ref_mean = dataForReferencing(2).data - dataForReferencing(2).mean;

dataForReferencing(3).data = (data_temp(gridCoor_3,:));
dataForReferencing(3).mean = nanmean(dataForReferencing(3).data,1);
dataForReferencing(3).data_ref_mean = dataForReferencing(3).data - dataForReferencing(3).mean;

dataForReferencing(4).data = (data_temp(gridCoor_4,:));
dataForReferencing(4).mean = nanmean(dataForReferencing(4).data,1);
dataForReferencing(4).data_ref_mean = dataForReferencing(4).data - dataForReferencing(4).mean;

data_ref_mean(gridCoor_1,:) = dataForReferencing(1).data_ref_mean;
data_ref_mean(gridCoor_2,:) = dataForReferencing(2).data_ref_mean;
data_ref_mean(gridCoor_3,:) = dataForReferencing(3).data_ref_mean;
data_ref_mean(gridCoor_4,:) = dataForReferencing(4).data_ref_mean;

%set the outlier electrodes to zero
data_ref_mean(outlierElectrodes,:) = 0;

end