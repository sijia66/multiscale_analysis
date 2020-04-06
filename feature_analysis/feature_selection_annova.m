
%% feature selection by annova

feature_p_values_all = ...
    fs_annova(s_z_trials_all_electrodes,trial_TargetAssigned);

%% examine one electrode

ei = 124;

feature_p_values = squeeze(feature_p_values_all(ei,:,:));
s_z_trials =  squeeze(s_z_trials_all_electrodes(ei,:,:,:));
ti = 1;
fi = 92;


figure
imagesc((feature_p_values'))
axis xy
colorbar
title('P Values')
xlabel('Time Points')
ylabel('Frequencies (Hz)')


feature_vec = squeeze(s_z_trials(:,ti,fi));
[p,~,stats] = anova1(feature_vec,target_vec);
xlabel('Target Direction')
ylabel('Averaged Z-score')

%% examine averages of electrodes
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all,3),2));

figure
imagesc(squeeze(mean(feature_p_values_all,2))')
ylabel('Frequencies')
xlabel('Electrodes')
axis xy
colorbar

figure
imagesc(squeeze(mean(feature_p_values_all,3))')
ylabel('Time Points')
xlabel('Electrodes')
axis xy
colorbar

figure
scatter(1:length(feature_p_values_mean),feature_p_values_mean)
xlabel('Electrodes')
ylabel('P-Values')
title('Comparing Behavior Correlation between Electrodes')

figure
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        squeeze(mean(mean(feature_p_values_all,3),2)));
imagesc(map')
scatterToMatrix_plot_label
title('P-values for ECOG Electrodes')

%% examine band specfic averages

plot_mean_freqBand(-log10(feature_p_values_all),freq,time_range_index)
xlabel('ECOG Electrodes')
ylabel('Log Transformed P-value')
set(gca, 'fontSize',22)


