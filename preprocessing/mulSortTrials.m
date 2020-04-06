function dataTarget = mulSortTrials(data,trial_TargetAssigned,Trials)
% This function sorts the data by their trials
% Inputs
%   data: trialN by Electrode by Time
%   trial_TargetAssigned: a column vector of size N*1 
% Output
%   dataTarget: a structure with Target
%

if size(data,1) ~= length(trial_TargetAssigned)
    error('The data structure length (%d) does not match the trial length (%d)'...
        ,size(data,1),length(trial_TargetAssigned))
end

%All checks out
disp('Start sorting!')

%initialize the matrix
dataTarget.target1 = [];dataTarget.target1_trials = Trials(trial_TargetAssigned == 1);
dataTarget.target2 = [];dataTarget.target2_trials = Trials(trial_TargetAssigned == 2);
dataTarget.target3 = [];dataTarget.target3_trials = Trials(trial_TargetAssigned == 3);
dataTarget.target4 = [];dataTarget.target4_trials = Trials(trial_TargetAssigned == 4);
dataTarget.target5 = [];dataTarget.target5_trials = Trials(trial_TargetAssigned == 5);
dataTarget.target6 = [];dataTarget.target6_trials = Trials(trial_TargetAssigned == 6);
dataTarget.target7 = [];dataTarget.target7_trials = Trials(trial_TargetAssigned == 7);
dataTarget.target8 = [];dataTarget.target8_trials = struct([]);

for ti = 1:size(data,1)
    switch uint8(trial_TargetAssigned(ti))
        case 1
            currLength = size(dataTarget.target1,1);
            dataTarget.target1(currLength+1,:,:) = data(ti,:,:);
        case 2
            currLength = size(dataTarget.target2,1);
            dataTarget.target2(currLength+1,:,:) = data(ti,:,:);
        case 3
            currLength = size(dataTarget.target3,1);
            dataTarget.target3(currLength+1,:,:) = data(ti,:,:);
        case 4
            currLength = size(dataTarget.target4,1);
            dataTarget.target4(currLength+1,:,:) = data(ti,:,:);
        case 5
            currLength = size(dataTarget.target5,1);
            dataTarget.target5(currLength+1,:,:) = data(ti,:,:);
        case 6
            currLength = size(dataTarget.target6,1);
            dataTarget.target6(currLength+1,:,:) = data(ti,:,:);
        case 7
            currLength = size(dataTarget.target7,1);
            dataTarget.target7(currLength+1,:,:) = data(ti,:,:);
        case 8
            currLength = size(dataTarget.target8,1);
            dataTarget.target8(currLength+1,:,:) = data(ti,:,:);
        otherwise
            disp('Something Wrong with this trial')
            disp(ti)
    end
end
disp('Done sorting the targets')

end