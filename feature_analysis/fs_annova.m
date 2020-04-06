function feature_p_values = fs_annova(s_z_trials,target_vec)
%This function calculates the p-values in using annova in selecting
%features
%Inputs:
%   s_z_trials_all_electrodes is the feature vector of the size:
%          by number of times by  number of frequencies by number of trials
%         
%   trial_TargetAssigned: target vector 
%           number of trials by 1
%Output
%   feature_p_values_all: un-corrected p-values of the size 
%                         number of electrodes by number of times by number
%                         of frequencies
% 
%Required Matlab packages
%Matlab statiscal learning toolboxes



[n_times, n_frequencies, n_trials] = size(s_z_trials);

%check if the number of trials of the feature matrix matches that of the
%target matrix
if length(target_vec) ~= n_trials
    error('Houston: Feature matrix'' number of trials does not match of the target vector')
end
%n_frequencies = 10; % just don


feature_p_values = zeros(n_times, n_frequencies);


for ti = 1:n_times
    for fi = 1:n_frequencies %only do lower diagnal eletrode pairs
        feature_vec = squeeze(s_z_trials(ti,fi,:));
        [p,~,stats] = anova1(feature_vec,target_vec, 'off');
        
        feature_p_values(ti, fi) = p;
        
    end
    if mod(ti,10) == 0
        %fprintf('Finished %d out of %d \n',ti,n_times)
    end
end


%disp('Finished Feature Selection by Annova')


end

