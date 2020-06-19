function [trial_order_ind, trial_order_val] = ...
                    sort_trial_by_target(trial_TargetAssigned)
% this function calculates the indices by which either trial or data
% it assumes the targets are 1,2,3,....N
% can be grouped together
% input:
%   trial_TargetAssigned: n by 1 vector with target directions
% output:
%   trial_order_ind: n by 1 vector where targets are grouped in increasing
%   order
%   trial_order_val: save as above with the actual targets in order
%example usage
%
                
MAX_TARGET_IND = max(trial_TargetAssigned);
trial_order_ind = [];
trial_order_val = [];
for i = 1:MAX_TARGET_IND
    temp = find(trial_TargetAssigned == i);
    trial_order_ind = [trial_order_ind; temp];
    trial_order_val = [trial_order_val; (double(i) * ones(size(temp)))];
end


end