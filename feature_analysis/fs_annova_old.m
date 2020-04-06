function feature_p_values_all = fs_annova_old(s_z_trials_all_electrodes,trial_TargetAssigned)
%This function calculates the p-values in using annova in selecting
%features
%Inputs:
%   s_z_trials_all_electrodes is the feature vector of the size
%          number of electrodes by number of trials by number of times by
%          number of frequencies 
%   trial_TargetAssigned: target vector 
%           number of trials by 1
%Output
%   feature_p_values_all: un-corrected p-values of the size 
%                         number of electrodes by number of times by number
%                         of frequencies
% 
%Required Matlab packages
%Matlab statiscal learning toolboxes



target_vec = trial_TargetAssigned;

[n_e, n_trials,  n_times, n_frequencies] = size(s_z_trials_all_electrodes);

%n_frequencies = 10; % just don 
feature_p_values = zeros(n_times, n_frequencies);

feature_p_values_all = [];

for ei = 1:n_e
    
    s_z_trials =  squeeze(s_z_trials_all_electrodes(ei,:,:,:));
    
    for ti = 1:n_times
        for fi = 1:n_frequencies
            feature_vec = squeeze(s_z_trials(:,ti,fi));
            [p,~,stats] = anova1(feature_vec,target_vec, 'off');
            
            feature_p_values(ti, fi) = p;
            
        end
    end
    
    feature_p_values_all(ei,:,:) = feature_p_values;
    if mod(ei,10) == 0
     fprintf('Finished electrode %d\n',ei)
    end
end

disp('Finished Feature Selection by Annova')


end

