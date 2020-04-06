function  add_t_f_axes(fig,t,f,trialInfo,num_ticks)
fig;
%add time points
t_length = length(t);
f_length = length(f);

%add movement start
xline(round(t_length * trialInfo.timeReachStart / (trialInfo.timeReachStart + trialInfo.tAfter)), 'LineWidth',2);

s_timePoints = t;
freq = f;

tick_step = round(t_length / num_ticks);
xticks(1:tick_step:t_length)
xticklabels(s_timePoints(1:tick_step:end))

tick_step = round(f_length / num_ticks);
yticks(1:tick_step:f_length)
yticklabels(round(freq(1:tick_step:f_length)))