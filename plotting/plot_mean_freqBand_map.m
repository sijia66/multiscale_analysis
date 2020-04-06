function plot_mean_freqBand_map(feature_p_values_all,P_chamber,freq,time_range_index)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% dependancies
%   cal_index _freq
%   scatterToMatrix
color_range = [0 5];


figure

freqRange = [0 4];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));

feature_p_values_mean(isnan(feature_p_values_mean)) = 1;
feature_p_values_mean =  -log10(feature_p_values_mean);


subplot(2,3,1)
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_mean);
imagesc(map')
%scatterToMatrix_plot_label
title('Theta (0 - 4 Hz)')
caxis(color_range)
axis xy


freqRange = [4 12];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
feature_p_values_mean(isnan(feature_p_values_mean)) = 1;
feature_p_values_mean =  -log10(feature_p_values_mean);

subplot(2,3,2)
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_mean);
imagesc(map')
%scatterToMatrix_plot_label
title('Delta (4 - 12 Hz)')
caxis(color_range)
axis xy

freqRange = [12 30];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
feature_p_values_mean(isnan(feature_p_values_mean)) = 1;
feature_p_values_mean =  -log10(feature_p_values_mean);

subplot(2,3,3)
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_mean);
imagesc(map')
%scatterToMatrix_plot_label
title('Beta (12 - 30 Hz)')
caxis(color_range)
axis xy

freqRange = [30 80];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
feature_p_values_mean(isnan(feature_p_values_mean)) = 1;
feature_p_values_mean =  -log10(feature_p_values_mean);

subplot(2,3,4)
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_mean);
imagesc(map')
%scatterToMatrix_plot_label
title('Gamma (30 - 80 Hz)')
caxis(color_range)
axis xy

freqRange = [80 300];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
feature_p_values_mean(isnan(feature_p_values_mean)) = 1;
feature_p_values_mean =  -log10(feature_p_values_mean);

subplot(2,3,5)
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_mean);
imagesc(map')
caxis(color_range)
%scatterToMatrix_plot_label
colorbar
title('High Gamma (80 - 300 Hz)')
axis xy
end

