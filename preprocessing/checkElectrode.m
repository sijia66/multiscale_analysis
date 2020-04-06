
%%
%query electrode properties

%assume each day, the depths of each electrode are the same
% what we ultimated want to do is to do a 3D rendering of the electrode
% position.
TrialN = 1;

%get the number of drives 
Ndrives = size(Trials(TrialN).Depth,2);
hold on

for di = 1:Ndrives
    
    P_chamber = drmap_coordTransform_microdrive2chamber(experiment, di);
    P_chamber = P_chamber'; % in mm
    
    acList = [experiment.hardware.microdrive(di).electrodes.acquisitionid];
    
    depthProfile = (Trials(1).Depth{1,di})'; % in micron
    depthProfileInMM = depthProfile / 1000;
    depthProfile_shortend = depthProfileInMM(acList)*(-1);
    
    scatter3(P_chamber(:,1),P_chamber(:,2), depthProfile_shortend)
end

hold off


%%
% here we will check the rms amplitude electrode
ti = 1;
data_temp = squeeze(data(ti,:,:));
data_temp_rms = rms(data_temp,2);

figure
scatter(1:length(data_temp_rms),data_temp_rms)
xlabel('Electrode Number')
ylabel('RMS Amplitude')
title(sprintf('RMS Amplitude for All Electrodes in Trial %d',ti))
 


