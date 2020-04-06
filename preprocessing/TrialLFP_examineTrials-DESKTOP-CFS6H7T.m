
global MONKEYDIR
MONKEYDIR = 'C:/Users/Orsborn Lab/OneDrive - UW/projects/Brain EEG/data';


%%
startToGo = [Trials.StartAq] - [Trials.Go];
goToOn = [Trials.Go] - [Trials.TargsOn];
onToStart = [Trials.TargsOn] - [Trials.ReachStart];
startToEnd = [Trials.ReachStart] - [Trials.ReachStop];

x = [startToGo';goToOn';onToStart';startToEnd'];
g = [ones(size(startToGo')); 2*ones(size(goToOn')); ...
    3*ones(size(onToStart'));4*ones(size(startToEnd'))];
boxplot(x,g,'Labels',{'StartAq to Go','Go to TargsOn','TargsOn to ReachStart',...
    'ReachStart to ReachEnd'})
ylabel('Time (ms)')


%%

global MONKEYDIR
MONKEYDIR = 'C:/Users/Orsborn Lab/OneDrive - UW/projects/Brain EEG/data';

trialNames = {'180321','180328','180329','180401','180408'...
                '180415'};

tNN = length(trialNames);

TrialN = 1;
di = 2;
magFac = 50; %bubble size
textShift = 0.2;

P_chamber_SC32 = ...
    drmap_coordTransform_microdrive2chamber(experiment, 'LM1_SC32_1');
P_chamber = P_chamber_SC32';

for tNi = 1:length(trialNames)
    tnTemp = trialNames{tNi};
    
    %get the trial structur
    trialDirectory = strcat(MONKEYDIR,'/',tnTemp,'/mat');
    cd(trialDirectory)
    
    load Trials.mat
    
    %get the number of drives
    Ndrives = size(Trials(TrialN).Depth,2);
    depthProfile = (Trials(1).Depth{1,di})'; % in micron
    
    depthProfiles(:,tNi) = depthProfile;
    
    disp(tnTemp)
end
% need to look

for tNi = 1:length(trialNames)
    
    depthProfile = depthProfiles(:,tNi);
    
    subplot(tNN,1,tNi)
    scatter(1:length(depthProfile),depthProfile)
    ylabel('Depth (um)')
    title(trialNames(tNi))
end
xlabel('SC32 Electrodes')


figure
% can plot this on a map

for tNi = 1:tNN
    
    depthProfile = depthProfiles(:,tNi);
    
    subplot(2,tNN/2,tNi)
    hold on
    for ei = 1:size(depthProfile,1)

        if depthProfile(ei) == 0
            scatter(P_chamber(ei,1),P_chamber(ei,2)...
                ,magFac, ...
                'MarkerFaceColor','g')
        elseif depthProfile(ei) > 0
            scatter(P_chamber(ei,1),P_chamber(ei,2)...
                ,magFac, ...
                'MarkerFaceColor','b')
        else
            scatter(P_chamber(ei,1),P_chamber(ei,2)...
                ,magFac, ...
                'MarkerFaceColor','y')
            
        end
        
        text(P_chamber(ei,1)+textShift,P_chamber(ei,2),num2str(ei))
    end
    hold off
    axis equal
    title(trialNames(tNi))
end