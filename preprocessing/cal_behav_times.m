function behav_struct = cal_behav_times(Trials)
%{
this function calculates the behaviour times relative to movement start and
assumes the sequence in behav_struct.labels
Inputs:
    Trials: num_trials assume the times are in ms

Outputs:
    behav_struct: a struct with the following field
    behav_struct.times relative to movement start in seconds
    behav_struct.times_std
    behav_struct.labels: hard coded labels for the sequence:
                        'Start acq' 'Go' 'Targ On' 'Reach Start' 'Reach Stop'

exmple usage

%}

startAcq = [Trials.StartAq] - [Trials.ReachStart];
targsOn = [Trials.TargsOn] - [Trials.ReachStart];
go = [Trials.Go] - [Trials.ReachStart];
reachStop = [Trials.ReachStop] - [Trials.ReachStart];

behav_times = [startAcq' go' targsOn' zeros(length(Trials),1) reachStop'];

behav_times_mean = nanmean(behav_times) / 1000; %convert to in s
behav_times_std = nanstd(behav_times) / 1000;   %convert to in s

behav_struct.times = behav_times_mean;
behav_struct.times_std = behav_times_std;
behav_struct.labels = {'Start acq' 'Go' 'Targ On' 'Reach Start' 'Reach Stop'};
