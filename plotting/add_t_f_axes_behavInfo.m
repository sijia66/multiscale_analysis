function  add_t_f_axes_behavInfo(fig,t,f,behav_times_mean,trialInfo,num_ticks)
%inputs
% behav_times: in s
% behav_times: 'start acq' 'go' 'target on' 'reach start' 'reach end'

fig;
%add time points
t_length = length(t);
f_length = length(f);

%add movement markers
delta_index = t_length / (t(end) - t(1)); %pixel for second
time_ReachStart_indx = abs(behav_times_mean(1) * delta_index);
behav_times_pixel = round(time_ReachStart_indx...
                    + delta_index*behav_times_mean);

%dirty fix of the begining and end of the display index array
behav_times_pixel(1) = 1;
behav_times_pixel(end) = length(t);
                
for i =  1:length(behav_times_pixel)
    xline(behav_times_pixel(i),'LineWidth',2);
end

%xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);

%add behav times

%calculate where to 



freq = f;

% tick_step = round(t_length / num_ticks);
% xticks(1:tick_step:t_length)
% xticklabels(s_timePoints(1:tick_step:end))


xticks(behav_times_pixel)
xticklabels(round(behav_times_mean,2))


tick_step = round(f_length / num_ticks);
yticks(1:tick_step:f_length)
yticklabels(round(freq(1:tick_step:f_length)))