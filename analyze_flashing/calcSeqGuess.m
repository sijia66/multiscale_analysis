function [Pos, TargEst, Time, TimeOff, trID, segID] = calcSeqGuess(Trials)

%[Pos, TargEst, Time, TimeOff, trID, segID] = calcSeqGuess(Trials)
%
%function to compile info (position, cloests target, times) of guesses (touches prior to the correct touch) in each segment of a sequene trial.
%
%inputs: Trials (struct, formatted according to Trials framework. Must include field '.seqEventsParsed' created by procSequence)
% 
%outputs:
%        Pos - (#guesses x 2) position of each guess touch in screen coordinates
%        TargEst - (#guesses x 2) position of the closest target in the array to the given touch position
%        Time    - (#guesses x 1) start time of touch, *relative to the
%                  go-cue for that given reach segment (or the startOn for
%                  the first segment which has no go-cue)*
%        TimeOff - (#guesses x 1) end time of touch. Same as 'Time',
%                  calculated relative to the go-cue for the given segment.
%        trID    - id for which trial each touch belongs to
%        segID   - id for which segment within the sequence each touch belongs to. 
%        
%
%A. Orsborn, 180118


%get only sequence learning trials
seqTr = find( [Trials.SequenceLearning]);

SEQ = [Trials(seqTr).seqEventsParsed];

TOUCH = [SEQ.touch];
%ON  = [SEQ.on];
GO  = [SEQ.go];


Pos        = [];
TargEst    = [];
Time       = [];
TimeOff    = [];
segID      = [];
trID       = [];
for i=1:length(seqTr)
    
    p    = [TOUCH(i).pos]; %positions of touches
    targ = [TOUCH(i).targEst]; 
    t    = [TOUCH(i).time];
    to   = [TOUCH(i).offTime];
    ids  = [TOUCH(i).segid];  %ID of segment
    
    if ~isempty(ids)
        
        %for each touch time, calculate it relative to the GO for that
        %segment. For first segment (which doesn't have a go, use startOn)
        refTimes = [Trials(seqTr(i)).StartOn [GO(i).time]];
        rerefmat = refTimes(ids+1); %segment ids are indexed from zero
        
        Pos     = cat(1, Pos, p);
        TargEst = cat(1, TargEst, targ);
        Time    = cat(1, Time, (t-rerefmat)');
        TimeOff = cat(1, TimeOff, (to-rerefmat)');
        segID   = cat(1, segID, ids);
        trID    = cat(1, trID, seqTr(i)*ones(size(ids)));
    end
    
end