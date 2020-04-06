
function dataTarget_matrix_min = mulSortTrials_matrix_min(dataTarget)
% This function finds the minimum trials of all the targets and outputs the
% same number of trials for each reach direction
% Inputs:
%   dataTarget
% Output
%   dataTarget_matrix_min


trialCount = [length(dataTarget.target1_trials);...
              length(dataTarget.target2_trials);...
              length(dataTarget.target3_trials);...
              length(dataTarget.target4_trials);...
              length(dataTarget.target5_trials);...
              length(dataTarget.target6_trials);...
              length(dataTarget.target7_trials);...
             ];
[trial_n,trial_n_target] = min(trialCount);
disp(sprintf('The minimum number of trials is %d with target %d', trial_n,trial_n_target));




disp('Assume there are only 7 recorded targets')
disp(sprintf('Saved the first %d trials for each target',trial_n))
for targi = 1:7
    eval(sprintf('dataTarget_matrix_min(targi,:,:,:) = dataTarget.target%d(1:trial_n,:,:);',targi))
end


end

