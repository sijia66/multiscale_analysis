function plot_trials(Trials)
% examine trials
startToGo = [Trials.Go] -[Trials.StartAq];
goToOn = [Trials.TargsOn] - [Trials.Go];
onToStart = [Trials.ReachStart] - [Trials.TargsOn];
startToEnd = [Trials.ReachStop] - [Trials.ReachStart];

x = [startToGo';goToOn';onToStart';startToEnd'];
g = [ones(size(startToGo')); 2*ones(size(goToOn')); ...
    3*ones(size(onToStart'));4*ones(size(startToEnd'))];
boxplot(x,g,'Labels',{'StartAq to Go','Go to TargsOn','TargsOn to ReachStart',...
    'ReachStart to ReachEnd'})
ylabel('Time (ms)')
set(gca,'fontSize', 6)
end