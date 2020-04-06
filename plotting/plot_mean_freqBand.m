function plot_mean_freqBand(feature_p_values_all,freq,time_range_index)
%This function averages frequency power in several bands 
%   Detailed explanation goes here


n_electrodes = size(feature_p_values_all);

figure


freqRange = [0 4];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:n_electrodes,feature_p_values_mean)

hold on
freqRange = [4 12];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:n_electrodes,feature_p_values_mean)


freqRange = [12 30];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:n_electrodes,feature_p_values_mean)

freqRange = [30 80];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:n_electrodes,feature_p_values_mean)

freqRange = [80 300];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:n_electrodes,feature_p_values_mean)

hold off
legend('Delta (0 - 4 Hz)',...
       'Alpha Theta (4 - 12 Hz)',...
       'Beta (12 - 30 Hz)',...
       'Gamma (30 - 80 Hz)',...
       'High Gamma (80 - 300 Hz)'...
        )

end

