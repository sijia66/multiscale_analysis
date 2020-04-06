function [Pos, targids] = calcSeqReachTarg(Trials)

%[Pos, targids] = calcSeqReachTarg(Trials)
%
%function to calculate the target position for each reach within the sequence task.
%
%inputs: Trials (struct, formatted according to Trials framework. Must include field '.seqEventsParsed' created by procSequence)
% 
%outputs:
%        Pos - (# trials x #targets x 2) matrix with hand position for each touch 
%                for each target within a sequence trial. Since sequence length may
%                vary across trials, #targets = max #targets across trials. Similarly, these calculations are not defined for non-sequenceLearning trials. 
%                Any missing data will be nans.
%        targids - (#trials x #targets) matrix with targid for each target
%                   within the sequence. center target = 0, 1st reach w/in sequence = 1, etc. See procSequence
%                   As with 'times' matrix, missing data will be nans. 
%
%A. Orsborn, 180118

nTr = size(Trials,2);

%get only sequence learning trials
seqTr = find( [Trials.SequenceLearning]);

SEQ = [Trials(seqTr).seqEventsParsed];

len = [SEQ.length];
ACQ = [SEQ.acq];
%ON  = [SEQ.on];

maxLen = max(len);

Pos = nan(nTr, maxLen*2, 2);
targids = nan(nTr, maxLen*2);
for i=1:length(seqTr)
    
    pos = [ACQ(i).targpos]; %position of target for each aquired target
    ids = [ACQ(i).targid];  %ID of each target being acquired (note that IDs are relative to the sequence, not the target array)

    
    %store in full matrix (dealing with variable sequence lengths)
    Pos(seqTr(i), 1:size(pos,1),:) = pos;
    targids(seqTr(i), 1:size(pos,1)) = ids;
end