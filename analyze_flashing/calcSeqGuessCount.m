function [count] = calcSeqGuessCount(Trials)

%[count] = calcSeqGuesCount(Trials)
%
%function to calculate the # of guesses (touches prior to the correct touch) in each segment of a sequene trial.
%
%inputs: Trials (struct, formatted according to Trials framework. Must include field '.seqEventsParsed' created by procSequence)
% 
%outputs:
%        counts - (# trials x #segments) matrix with guess count (# of touches)
%                for each segment within a sequence trial. #segments = 2*length of sequence (i.e. segment 1 = touches prior to acquiring initial center. segment 2 = touches prior to acquiring target1, etc.)
%                Since sequence length may vary across trials, #targets = max #targets across trials. 
%                Similarly, these calculations are not defined for non-sequenceLearning trials. 
%                Any missing data will be nans.
%
%A. Orsborn, 180118

nTr = size(Trials,2);

%get only sequence learning trials
seqTr = find( [Trials.SequenceLearning]);

SEQ = [Trials(seqTr).seqEventsParsed];

len = [SEQ.length];
TOUCH = [SEQ.touch];
%ON  = [SEQ.on];
%GO  = [SEQ.go];

maxLen = max(len);

count = nan(nTr, maxLen*2); % # segments in the sequence = 2*seqLen b/c of center cues.
for i=1:length(seqTr)
    
    %touchOn = [TOUCH(i).time]; %time when touched
    ids     = [TOUCH(i).segid];  %ID of segment
    
    if ~isempty(ids)
        ids2 = repmat(ids, 1, len(i)*2); %ids x #segments matrix
        segs = repmat(0:(len(i)*2-1), size(ids,1), 1); %segid matrix ([0 1 2 3 4...; 0 1 2 3 4...;])
        
        count(seqTr(i),1:len(i)*2) = sum(ids2==segs, 1);
    else
        count(seqTr(i),1:len(i)*2) = zeros(1, len(i)*2);
    end
end