function [TargPos, TargColor, CenterPos, CenterColor] = calcSeqTargArray(Trials)

%[Pos, targids] = calcSeqTargArray(Trials)
%
%function to calculate the target array properties for sequence trials.
%
%inputs: Trials (struct, formatted according to Trials framework. Must include field '.SEQ' created by procEvents)
% 
%outputs:
%        TargPos - (# trials x #targets x 2) matrix with target array position for each trial 
%                for each target within a sequence trial. Since sequence length may
%                vary across trials, #targets = max #targets across trials. Similarly, these calculations are not defined for non-sequenceLearning trials. 
%                Any missing data will be nans.
%        TargColor - (#trials x #targets x 3) matrix with target color values (RGB) for each trial, each targettargid for each target
%                   within the sequence. As with TargPos, missing data will be nans. 
%        CenterPos - (#trials x 2) matrix with center target location for each trial
%                     missing values are nan
%        CenterColor - (#trials x 3) matrix with center color values for each trial. 
%                     missing values are nan
%
%A. Orsborn, 180118

nTr = size(Trials,2);

%get only sequence learning trials
seqTr = find( [Trials.SequenceLearning]);

SEQ = [Trials(seqTr).SEQ];

setsz = [SEQ.setSize];
%ACQ = [SEQ.acq];
%ON  = [SEQ.on];
%GO  = [SEQ.go];

maxLen = max(setsz);

TargPos     = nan(nTr, maxLen, 2);
TargColor   = nan(nTr, maxLen,3);
CenterPos   = nan(nTr, 2);
CenterColor = nan(nTr, 3);

for i=1:length(seqTr)
    
    tpos = [SEQ(i).targLocations]; %positions of target arrays
    tcol = [SEQ(i).targColors]; %colors of target arrays
    cpos = [SEQ(i).centerLocations]; %center pos
    ccol = [SEQ(i).centerColors];    %center color
    
    %store in full matrix (dealing with variable sequence lengths)
    TargPos(seqTr(i), 1:size(tpos,2),:) = tpos';
    TargColor(seqTr(i), 1:size(tcol,2),:) = tcol';
    
    CenterPos(seqTr(i),:) = cpos;
    CenterColor(seqTr(i),:) = ccol;
end

