%% try out the example in matlab

load hogg
hogg

[p,~,~] = anova1(hogg,'off');

%% try out our own example

strength = [82 86 79 83 84 85 86 87 74 82 ...
            78 75 76 77 79 79 77 78 82 79];
alloy = {'st','st','st','st','st','st','st','st',...
         'al1','al1','al1','al1','al1','al1',...
         'al2','al2','al2','al2','al2','al2'};
     
     [p,~,stats] = anova1(strength,alloy,'off');
     
     
%%  a N by 1 vector and a target reach veotor

ti =  1;
fi = 1;

feature_vec = squeeze(s_z_trials(:,ti,fi));
target_vec = trial_TargetAssigned;

[p,~,stats] = anova1(feature_vec,target_vec);

% from this plot, we are separating some groups, it might be useful 
% to write our own analysis code
% what are other ways to carry out this  analysis

%% after this practice, we can write a double for loop and figure out the 
% the features
% we can also do electrodes ! 
target_vec = trial_TargetAssigned;

[n_trials, n_times, n_frequencies] = size(s_z_trials);

%n_frequencies = 10; % just don 
feature_p_values = zeros(n_times, n_frequencies);

for ti = 1:n_times
    for fi = 1:n_frequencies
        feature_vec = squeeze(s_z_trials(:,ti,fi));
        [p,~,stats] = anova1(feature_vec,target_vec, 'off');
        
        feature_p_values(ti, fi) = p;
       
    end
    disp(ti)
end

disp('Done!')

figure
imagesc((feature_p_values'))
axis xy
colorbar
title('P Values')
xlabel('Time Points')
ylabel('Frequencies ')



%% we can look several features if they make sense
ti = 3;
fi = 21;
feature_vec = squeeze(s_z_trials(:,ti,fi));
[p,~,stats] = anova1(feature_vec,target_vec);

ti = 1;
fi = 92;
feature_vec = squeeze(s_z_trials(:,ti,fi));
[p,~,stats] = anova1(feature_vec,target_vec);
%% next thing we can do is to compare electrodes 
% 


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

disp('Done!')
%% look at the p-values for these electrodes
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
scatter(1:211,feature_p_values_mean)
xlabel('Electrodes')
ylabel('P-Values')
title('Comparing Behavior Correlation between Electrodes')

figure
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        squeeze(mean(mean(feature_p_values_all,3),2)));
imagesc(map')
scatterToMatrix_plot_label
title('P-values for ECOG Electrodes')


%% look at a single electrode

ei = 124;

feature_p_values = squeeze(feature_p_values_all(ei,:,:));

figure
imagesc((feature_p_values'))
axis xy
colorbar
title('P Values')
xlabel('Time Points')
ylabel('Frequencies (Hz)')

s_z_trials =  squeeze(s_z_trials_all_electrodes(ei,:,:,:));
ti = 1;
fi = 92;
feature_vec = squeeze(s_z_trials(:,ti,fi));
[p,~,stats] = anova1(feature_vec,target_vec);
xlabel('Target Direction')
ylabel('Averaged Z-score')

%% compile into one function 
ei = 1;

s_z_trials_one_electrode = (s_z_trials_all_electrodes(ei,:,:,:));
feature_p_values_all_temp  = ...
    fs_annova(s_z_trials_one_electrode,trial_TargetAssigned);

imagesc(squeeze(feature_p_values_all_temp))


%% along the same vein, we are gonna look at smallest p-values for each electrode

feature_p_values_min = squeeze(min(min(feature_p_values_all,[],3),[],2));

figure
scatter(1:211,feature_p_values_min)
xlabel('ECOG Electrodes')
ylabel('P-Values')
title('Comparing Behavior Correlation between Electrodes (Min P-value)')

figure
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_min);
imagesc(map')
scatterToMatrix_plot_label
title('Min P-values for Each ECOG Electrode')

%% look at band specific features 
%what are the bands to look for 
time_range_index = 3;

freqRange = [0 4];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));


figure
plot(1:211,feature_p_values_mean)

hold on
freqRange = [4 12];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:211,feature_p_values_mean)


freqRange = [12 30];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:211,feature_p_values_mean)

freqRange = [30 80];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:211,feature_p_values_mean)

freqRange = [80 300];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));
plot(1:211,feature_p_values_mean)

hold off
legend('Delta (0 - 4 Hz)',...
       'Alpha Theta (4 - 12 Hz)',...
       'Beta (12 - 30 Hz)',...
       'Gamma (30 - 80 Hz)',...
       'High Gamma (80 - 300 Hz)'...
        )
xlabel('ECOG Electrodes')
ylabel('Averaged P-value')
set(gca, 'fontSize',22)
%% we are going to convert this to 
time_range_index = 3;


color_range = [0 0.5];

freqRange = [0 4];
freq_range_index  = cal_index_freq(freq ,freqRange(1),freqRange(2));
feature_p_values_mean = squeeze(mean(mean(feature_p_values_all(:,time_range_index,freq_range_index ),3),2));


figure
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

subplot(2,3,5)
[map,mapInfo] = scatterToMatrix(P_chamber, ...
                        feature_p_values_mean);
imagesc(map')
caxis(color_range)
%scatterToMatrix_plot_label
colorbar
title('High Gamma (80 - 300 Hz)')
axis xy

%% compile into one function

plot_mean_freqBand(-log10(feature_p_values_all),freq,time_range_index)
xlabel('ECOG Electrodes')
ylabel('Log Transformed P-value')
set(gca, 'fontSize',22)

%% plot the same thing, but do a log transformation 

time_range_index = 3;
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
%% compile the spatial mapping into a function
plot_mean_freqBand_map(feature_p_values_all,P_chamber_ECOG_3,freq,time_range_index)


%% analyzing for SC32 electrodes

feature_p_values_all_SC32 = ...
    fs_annova(s_z_trials_all_electrodes_SC32,trial_TargetAssigned);


%% analyze one electrode
ei = 6;

feature_p_values = squeeze(feature_p_values_all_SC32(ei,:,:));
s_z_trials =  squeeze(s_z_trials_all_electrodes(ei,:,:,:));
ti = 3;
fi = 20;


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

%% look at averages of electrodes for SC32
% do a little transformation here
depthProfile = (Trials(1).Depth{1,drive_idx})'; % in micron
loweredElectrodes  = find(depthProfile > 0);
unlowered_electrodes =  find(depthProfile == 0);

feature_p_values_all_SC32(unlowered_electrodes,:,:) = 1;

plot_mean_freqBand(-log10(feature_p_values_all_SC32),freq,time_range_index)
xlabel('Lowered SC32 Electrodes')
ylabel('Log Transformed P-value')
set(gca, 'fontSize',22)

%% plot the SC32 electrodes 

plot_mean_freqBand_map(feature_p_values_all_SC32,P_chamber_SC32,freq,time_range_index)
title('SC32')

