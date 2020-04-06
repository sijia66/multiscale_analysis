function [trLfpData,trialInfo] = load_data(Trials, experiment,trialInfo,MONKEYDIR)
%this function loads combined ECOG and LFP
%   outputs: trials * number of electrodes * time
%       


    %first load ECOG
    driveNameAnalyze = 'LM1_ECOG_3';
    drive_idx = 1; %1 for ECOG and 2 for SC3
    dataParam.eList = [experiment.hardware.microdrive(drive_idx).electrodes.channelid];
    data_ECOG = trialLfp(Trials,driveNameAnalyze,dataParam.eList,[],trialInfo.trig,[trialInfo.tBefore,trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);
    trialInfo.ECOG_indices = 1:size(data_ECOG,2);
    
    
    %then we load LFP
    driveNameAnalyze = 'LM1_SC32_1';
    drive_idx = 2; %1 for ECOG and 2 for SC32
    dataParam.eList = [experiment.hardware.microdrive(drive_idx).electrodes.channelid];
    data_SC32 = trialLfp(Trials,driveNameAnalyze,dataParam.eList,[],trialInfo.trig,[trialInfo.tBefore,trialInfo.tAfter],trialInfo.lfpType,MONKEYDIR);
    trialInfo.SC32_indices = (size(data_ECOG,2) + 1):(size(data_ECOG,2) + size(data_SC32,2));
    
    trLfpData = cat(2,data_ECOG, data_SC32);
    
    clear data_SC32
    clear data_ECOG

end

