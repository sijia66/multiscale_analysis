function TargPosEst = calcSeqClosestTarget(acqPos, acqTargID, TargPos, CenterPos)
% 
%TargPosEst = calcSeqClosestTarget(acqPos, TargPos, CenterPos)
%
%
%This function estimates the position of the presented target being acquired in the sequence task. 
%This is particularly useful as a work-around if matching of labview and event-based trials proves tricky. 
%
%inputs:acqPos - (# trials x # movements w/in sequence x 2) matrix of touch positions at time of target acquisition for a sequence
%                 Can be generated with calcSeqReachPos
%       TargPos - (# trials x # targets in array x 2) matrix of target array
%                  positions for each sequence trial 
%                  Can be generated with calcSeqTargArray
%       CenterPos - (# trials x 2) matrix of center positions for each sequence trial
%                   Can be generated with calcSeqTargArray
%
%outputs: TargPosEst - (# trials x #movements w/in sequence x 2) matrix
%                       with closest target position of ALL presented targets (across trials)  for
%                       each target acquisition
%
%note: does NOT assume that hand-scaling has been done. Normalizes touch
%positions and target positions to be zero-mean, unit variance before
%looking for position matches. 
%
%
%A. Orsborn, 181024
%        

%%% find presented target position for each touch
%This is a work-around for poor matches between labview trials and event-based trials...

%find unique targets within the arrays shown in session
TargPos = reshape(TargPos, [size(TargPos,1)*size(TargPos,2),2]);
UniqueTargPos = unique(TargPos, 'rows');
UniqueCPos    = unique(CenterPos, 'rows');
nCenters      = size(UniqueCPos,1);
UniqueTargPos = [UniqueCPos; UniqueTargPos];
nTarg = size(UniqueTargPos,1);

% %scale targ array to be zero-mean, range 1
% Mx = max(UniqueTargPos,[],1);
% Mn = min(UniqueTargPos, [],1);
% M = mean(UniqueTargPos, 1);
% UniqueTargPos_scaled = (UniqueTargPos - repmat(M, nTarg,1))./ repmat((Mx - Mn), nTarg,1);
% 
UniqueTargPos_scaled = UniqueTargPos;

%similarly get re-scaling for touch positions
tmp = reshape(acqPos, [size(acqPos,1)*size(acqPos,2), 2]); %re-shape to tr*touches x 2 
%remove any nans
Mx_tch = max(tmp,[],1);
Mn_tch = min(tmp, [],1);
M_tch = nanmean(tmp, 1);


%for each touch, find target position that's closest to touch.

[nTr, nTouch, ~] = size(acqPos);
TargPosEst = nan(size(acqPos));
for tr=1:nTr
    for to=1:nTouch
        
%         %scale curr position
%         T = (sq(acqPos(tr,to,:))' - M_tch)./ (Mx_tch - Mn_tch);
        T = sq(acqPos(tr,to,:))';
        
        if ~any(isnan(T))
        T = repmat(T,nTarg,1); %replicate # targets x 2
        
        %find position w/ closest distance
        dist = sum((UniqueTargPos_scaled-T).^2,2);
        [~,tInd] = min( dist);
        
        %clean up noisy touches that sometimes seem closer to the center
        %(mostly happens with sliding for upper targets)
        if acqTargID(tr,to)~=0 && tInd<=nCenters %not a center reach but found center is closest
            dist(tInd) = inf; %find next-closest target that isn't the center
            [~,tInd] = min( dist);
            disp('Warning. Force-fixing Target estimate assignments based on targID')
        elseif acqTargID(tr,to)==0 && tInd>nCenters %if it's a center reach but found a non-central target
            dist(nCenters+1) = inf; %force it to find a center target
            [~,tInd] = min( dist);
            disp('Warning. Force-fixing Target estimate assignments based on targID')
        end

        
        TargPosEst(tr,to,:) = UniqueTargPos(tInd,:);
        else
            TargPosEst(tr, to, :) = nan(1,1,2);
        end
    end
end
   


