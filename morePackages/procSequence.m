function Events = procSequence(day, rec, debugging)
%This function examines task events and touch events to parse events within
%sequence trials
%This includes: time of touches w/in a sequence (may be multiple)
%               closest target to each detected touch




global MONKEYDIR

touchAvgWin = 100;

if ~exist('debugging', 'var')
    debugging = 0;
end

if exist('rec', 'var') && ~isempty(rec)
    recs = parseRecInput(rec);
else
    recs = dayrecs(day, MONKEYDIR);
end

nRec = length(recs);

%loop across recordings
for iR=1:nRec
    
    fname_prefix = [MONKEYDIR '/' day '/' recs{iR} '/rec' recs{iR} ];
    
    %double-check that events file exists. If not, skip recording.
    if ~exist([fname_prefix, '.Events.mat'], 'file')
        disp(['Events file for rec' recs{iR} ' does not exist. Skipping.'])
        continue
    else
        
        %load events, hand calibration data, experiment def
        load([fname_prefix '.Events.mat']);
        if exist([fname_prefix '.HandScale.mat'], 'file')
            load([fname_prefix '.HandScale.mat']);
        else
            HandScale = [];
        end
        load([fname_prefix '.experiment.mat']);
        
        %get format of hand data file
        if isfield(experiment,'hardware')
            format = getFileFormat('Broker');
        elseif isfile(['rec' rec '.experiment.mat'])
            load(['rec' rec '.experiment.mat'])
            format = getFileFormat('Broker');
        elseif isfile(['rec' rec '.Rec.mat']);
            load(['rec' rec '.Rec.mat']);
            if ~isfield(Rec,'BinaryDataFormat')
                Rec.BinaryDataFormat = 'short';
            end
            format = Rec.BinaryDataFormat;
        end
        
        %clear out any old sequence procing
        if isfield(Events, 'seqEventsParsed')
            Events = rmfield(Events, 'seqEventsParsed');
        end
        
        
        %loop through each trial
        for iTr = 1:length(Events.Trial)
            clear seq
            
            %only look at sequence trials
            if Events.SequenceLearning(iTr)
                

                %load all events that occurred for given sequence, sequence
                %info
                evs      = Events.seqEvents{iTr};
                disp_evs = Events.seqDispEvents{iTr};
                targPos  = Events.SEQ(iTr).targLocations;
                cPos     = Events.SEQ(iTr).centerLocations;
                seqLen   = Events.SequenceLength(iTr); %SEQ.length seems to be off?
                
                %load hand data for this trial. align to first event in a sequence (startOn), load from [0 lastEvent])
                %note: startOn field has been adjusted via displayEvents.
                len   = evs(end,3) - evs(1,3) +1;
                Hand = loadhnd(fname_prefix, Events, iTr, 'StartOn', [0 len], 1, HandScale, format);
                
                trStartTime = Events.StartOn(iTr);
                
                
                %%%%find discrete touches
                %set thresholds
                [~, bn] = hist(Hand(1,:));
                xTouchThresh = bn(2);
                
                [~,bn] = hist(Hand(2,:));
                yTouchThresh = bn(2);
                
                %only consider touch when both x, y above threshold
                touchingScreen = (Hand(1,:) > xTouchThresh) & (Hand(2,:) > yTouchThresh);
                
                %get touch start times and end times
                touchStarts = find(diff(touchingScreen)==1);
                touchEnds   = find(diff(touchingScreen)==-1);
                
                %deal with edge-cases
                if ~isempty(touchStarts) && ~isempty(touchEnds)
                    %if touching screen at start of trial and missed start
                    if touchEnds(1)<touchStarts(1)
                        touchStarts = [1 touchStarts];  %#ok<AGROW>
                    end
                    if length(touchEnds)<length(touchStarts)
                        touchEnds(end+1) = len; %#ok<AGROW> %if still touching screen at time of reward, consider reward time as end of touch
                    end
                else
                    %if only one touch, do same thing
                    if isempty(touchStarts) && length(touchEnds)==1
                        touchStarts(1) = 1;
                    end
                    if isempty(touchEnds) && length(touchStarts)==1
                        touchEnds(1) = len;
                    end
                end
                
                
                %may want to add a 'blink' to ignore touches starting <10ms after prev one
                %ends (??)--alo 180930
                nTouch = length(touchStarts);
                
                %keep properties of each touch within the trial
                seqTouchStart = touchStarts + trStartTime;
                seqTouchEnd   = touchEnds + trStartTime;
                seqTouchPos   = zeros(2,nTouch);
                seqTouchTargEst   = zeros(2,nTouch);
                seqTouchTargEstID = zeros(1,nTouch);
                seqTouchTargErr   = zeros(2,nTouch);
                
                fullTargList = [cPos targPos]; %center + target array
                for t=1:nTouch
                    last = min(touchEnds(t)-touchStarts(t), touchAvgWin);
                    seqTouchPos(:,t) = median( Hand(:,(1:last)+touchStarts(t)),2); %Take average at beginning of touch to deal with sliding. Want the initial touch position. 
                    
                    %estimate closest target in array (useful for guesses)
                    d = fullTargList - repmat(seqTouchPos(:,t), 1, size(targPos,2)+size(cPos,2));
                    d2 = sqrt( sum(d.^2,1));
                    [~, closest] = min(d2);
                    
                    seqTouchTargEst(:,t) = fullTargList(:,closest);
                    seqTouchTargEstID(t) = closest-1; %center = 0, targets in list
                    seqTouchTargErr(:,t) = d(:, closest);
                    
                end %loop touches
                
                
                
                %%%%%%%
                %%now parse codes out across sequence trial
                onTimes  = nan(1, seqLen*2); %time of target coming on
                onID     = nan(1, seqLen*2); %ID of target for each on time
                acqTimes = nan(1, seqLen*2); %times of each target acquisition within sequence
                acqOffTimes = nan(1, seqLen*2); %times of touch ends (hand off screen) within sequence
                acqPos   = nan(seqLen*2, 2); %position of touch for target acquisition
                acqTargPos = nan(seqLen*2, 2); %position of target aiming for
                acqID    = nan(1, seqLen*2); %ID of target for each acq time
                goTimes  = nan(1, seqLen*2-1); %time of go-cue for each target
                goID     = nan(1, seqLen*2-1); %ID of target for go time
                touchTimes    = []; %time (start time) of each 'extra' touche outside of acquiring a target
                touchTimesEnd = []; %time (end time) of each 'extra' touch outside of accquiring a target
                touchPos      = []; %position of each touch
                touchTargEst  = []; %closest target to touch
                touchID       = []; %ID of trial interval each touch belongs to
                
                
                
                %%%first center target--center on, acquiring center, guesses in between
                onTimes(1) = Events.StartOn(iTr); %first center, array on
                onID(1)    = 0; %0 = center
                
                if any(evs(:,2)==2)
                    acqTimes(1) = evs( evs(:,2)==2, 3 ); %center acquired code
                    acqID(1)   = 0; %0=center
                    
                    %find touch near acqTime to log touch position for initial
                    %center touch
                    touchInds = seqTouchStart>evs(1,3) & seqTouchStart<=acqTimes(1)+100; %all touches w/in segment
                    correctTargTouch = find( seqTouchTargEstID==0 & touchInds, 1, 'last');
                    %                 a = seqTouchStart(correctTargTouch); %use this if want to
                    %                 use touch-based estimate of center acquisition rather
                    %                 than event-code
                    %                 if isempty(a)
                    %                     a = nan;
                    %                 end
                    if ~isempty(correctTargTouch)
                        acqPos(1,:)     = seqTouchPos(:,correctTargTouch);
                        acqTargPos(1,:) = seqTouchTargEst(:,correctTargTouch);
                        acqOffTimes(1)  = seqTouchEnd(correctTargTouch);
                    else
                        acqPos(1,:)     = nan(1,2);
                        acqTargPos(1,:) = nan(1,2);
                        acqOffTimes(1)  = nan;
                        
                    end
                    
                else
                    acqTimes(1) = nan;  %set nan for center
                    acqOffTimes(1) = nan;
                    acqID(1)    = 0;
                    acqPos(1,:) = nan(1,2);
                    acqTargPos(1,:) = nan(1,2);
                end
                
                
                
                %log any touches occurring during interval between trial start and center acquired
                buff = 20; %if w/in 10ms, ignore. Probably just jitter in touch recognition
                if ~isempty(acqTimes) %seq.c1Acq)
                    g = seqTouchStart < acqTimes(1) - buff; %seq.c1Acq - buff;
                    %if center not acquired, log all touches between start on and
                    %error
                else
                    g = seqTouchStart < evs(evs(:,2)==8, 3); %error code
                end
                
                if any(g) %just avoids empties w/ weird sizes
                    touchTimes = cat(2, touchTimes, seqTouchStart(g));
                    touchTimesEnd = cat(2, touchTimesEnd, seqTouchEnd(g));
                    touchPos = cat(1, touchPos, seqTouchPos(:,g)');
                    touchTargEst = cat(1, touchTargEst, seqTouchTargEst(:,g)');
                    touchID  = cat(1, touchID,  zeros(sum(g),1)); %all these touches are in segment 0 (pre-center)
                end
                
                %%%subsequent targets
                goInds = find(evs(:,2)==4); % # of go cues = #targets presented after initial center cue
                %need to loop through sequenceLength so that event codes are
                %created for entire sequence even if trial aborted early (e.g.
                %if got an error at target 1 on a length 2 sequence, still want
                %to generate nan events for targ2
                nReaches = seqLen*2-1;
                for t=1:nReaches;
                    
                    
                    if mod(t,2)==1
                        %targInd = (t+1)/2;      %target indexing is every other, starting w/ 1 (1st go = 1st targ; 3rd go = 2nd targ)
                        targID  = (t+1)/2;      %ID for target in list
                    else
                        %targInd = (t+2)/2;      %center indexing is every other, starting w/ 2 (2nd go = 2nd c; 4th go = 3rd c, etc.)
                        targID  = 0;            %ID for center is 0 (see coding of touch targ est above)
                    end
                    
                    %if ran out of reach segments, set ons and acqs to nan
                    if t>length(goInds)
                        
                        onTimes(t+1) = nan;    %note on and acq shifted by 1 b/c also exist for initial center targ
                        onID(t+1)   = targID;
                        
                        acqTimes(t+1) = nan;
                        acqOffTimes(t+1) = nan;
                        acqPos(t+1,:) = nan(1,2);
                        acqTargPos(t+1,:) = nan(1,2);
                        acqID(t+1)   = targID;
                        
                        goTimes(t) = nan;
                        goID(t)    = targID;
                        
                        %if still have data, parse reach events
                    else
                        
                        %sub-select events just within this reach segment (i.e. between go cues or go cue to end)
                        ind1 = goInds(t);
                        if t<length(goInds)
                            ind2 = goInds(t+1)-1;
                            ind3 = goInds(t+1);
                        else
                            ind2 = size(evs,1);
                            ind3 = ind2;
                        end
                        evs_segment = evs(ind1:ind2,:);
                        maxTsearch  = evs(ind3,3); %limit display event/touch searches w/in time range between current go cue and next go-cue
                        
                        
                        
                        %go cue time is next display event after go-cue code
                        g = disp_evs(find(disp_evs(:,2)>evs_segment(1,3) & disp_evs(:,2)>onTimes(1), 1, 'first'),2);
                        goTimes(t) = g;
                        goID(t)   = targID;
                        
                        
                        %if visual cue for target occurred for this target,
                        %there will be another display event after this one,
                        %but less than next go cue/end of trial
                        cueTime = disp_evs(find(disp_evs(:,2)>g & disp_evs(:,2)<maxTsearch, 1, 'first'),2);
                        if isempty(cueTime)
                            cueTime = maxTsearch; %set cue time for next target--this just assures that acq time - cuetime is negative, which is used for subsequent calc.
                        end
                        onTimes(t+1) = cueTime;
                        onID(t+1)   = targID;
                        
                        
                        %%%find time when touched the target
                        touchInds = seqTouchStart>evs_segment(1,3) & seqTouchStart<maxTsearch; %all touches w/in segment
                        
                        %find LAST touch of the current targetID--assume this
                        %is the acquisition of the target
                        correctTargTouch = find( seqTouchTargEstID==targID & touchInds, 1, 'last');
                        a = seqTouchStart(correctTargTouch);
                        if isempty(a)
                            %if can't find touch via target position, use
                            %touch closest to acquire code in segment
                            
                            %only do for non-error trials
                            if isempty(evs_segment(evs_segment(:,2)==8,3))
                                
                                acqCodeTime = evs_segment( evs_segment(:,2)==5, 3);
                                if isempty(acqCodeTime)
                                    acqCodeTime = evs_segment( evs_segment(:,2)==3,3);
                                elseif length(acqCodeTime)>1
                                    disp(['Warning. Event codes have unusual structure. Check:' day ', rec: ' recs{iR} ', trial: ' num2str(iTr)])
                                    acqCodeTime = acqCodeTime(1);
                                end
                                
                                %only search touches within range of go-cue
                                %and end of this segment
                                touchInds = find(seqTouchStart>=goTimes(t)  & seqTouchStart<= maxTsearch);
                                touches = seqTouchStart(touchInds);
                                %if NO touches start within this segment, probably have some 'sliding' going on. So
                                %use the acqCode timing rather than touch start time
                                if isempty(touches)
                                    a = acqCodeTime;
                                    correctTargTouch = nan; %don't use touch start positions.
                                    disp(['Sliding detected...' day ', rec: ' recs{iR} ', trial: ' num2str(iTr) ', seg: ' num2str(t)])
                                else
                                    
                                    [tmp, ind] = min( abs(touches - acqCodeTime) ); %min( abs(seqTouchStart - acqCodeTime) ) ;
                                    if tmp>200
                                        disp(['Warning. Target acquisition timing questionable on day:' day ', rec: ' recs{iR} ', trial: ' num2str(iTr) ', seg: ' num2str(t)])
                                    end
                                    a = touches(ind);
                                    correctTargTouch = touchInds(ind);
                                end
                            else
                                %if not a correct segment, don't log a touch
                                a = nan;
                            end
                        end
                        acqTimes(t+1)     = a;
                        acqID(t+1)        = targID;
                        if ~isempty(correctTargTouch)
                            if ~isnan(correctTargTouch)
                                acqPos(t+1,:)     = seqTouchPos(:,correctTargTouch);
                                acqTargPos(t+1,:) = seqTouchTargEst(:,correctTargTouch);
                                acqOffTimes(t+1)  = seqTouchEnd(correctTargTouch);
                            else
                                %if there was sliding, need to calc touch
                                %pos etc
                                acqPos(t+1,:) = mean(Hand(:,(-20:20)+acqTimes(t+1)-trStartTime), 2);
                                d = fullTargList - repmat(acqPos(t+1,:)', 1, size(fullTargList,2));
                                d2 = sqrt( sum(d.^2,1));
                                [~, closest] = min(d2);
                                
                                acqTargPos(t+1,:) = fullTargList(:,closest);
                                acqOffTimes(t+1)  = seqTouchEnd( find(seqTouchEnd>a, 1, 'first'));
                                
                            end
                                
                        else
                            acqPos(t+1,:) = nan(1,2);
                            acqTargPos(t+1,:) = nan(1,2);
                            acqOffTimes(t+1)  = nan;
                        end
                        
                        %log any touches occurring during interval between trial start and target acquired (or error)
                        if ~isnan(a)
                            endTime = a; %acquired target
                        else
                            endTime = evs_segment(evs_segment(:,2)==8, 3); %error code
                            if isempty(endTime)
                                endTime=evs_segment(end,3);
                            end
                        end
                        g = seqTouchStart > evs_segment(1,3) & seqTouchStart< endTime - buff;
                        
                        if any(g)
                            touchTimes = cat(2, touchTimes, seqTouchStart(g));
                            touchTimesEnd = cat(2, touchTimesEnd, seqTouchEnd(g));
                            touchPos = cat(1, touchPos, seqTouchPos(:,g)');
                            touchTargEst = cat(1, touchTargEst, seqTouchTargEst(:,g)');
                            touchID  = cat(1, touchID,  t*ones(sum(g),1)); %reach segments within sequence
                        end
                        
                    end %check for reach segments
                    
                end %loop through go-cues (reach segments)

                %put together into structure
                clear seq
                seq.acq.time      = acqTimes;
                seq.acq.offTime   = acqOffTimes;
                seq.acq.pos       = acqPos;
                seq.acq.targPos   = acqTargPos;
                seq.acq.targid    = acqID;
                seq.go.time       = goTimes;
                seq.go.targid     = goID;
                seq.on.time       = onTimes;
                seq.on.targid     = onID;
                seq.touch.time    = touchTimes;
                seq.touch.offTime = touchTimesEnd;
                seq.touch.pos     = touchPos;
                seq.touch.targEst = touchTargEst;
                seq.touch.segid   = touchID;
                
                %log reward on, error
                seq.rewardOn = evs(evs(:,2)==7,3);
                if isempty(seq.rewardOn)
                    seq.rewardOn = nan;
                end
                seq.errorOn  = evs(evs(:,2)==8,3);
                if isempty(seq.errorOn)
                    seq.errorOn = nan;
                end
                
                seq.length = Events.SequenceLength(iTr);
                seq.rand   = Events.SEQ(iTr).rand;
                seq.index  = Events.SEQ(iTr).index;
                
                Events.seqEventsParsed(iTr) = seq;
                
                
                %load display, plot things out (for debugging)
                %acqRelCue = acqTimes - onTimes
                    
                if debugging %& any(acqRelCue==0)
                    if ~exist('debugfig', 'var')
                        debugfig = figure;
                    else
                        figure(debugfig);
                        clf(debugfig)
                        hold off
                    end
                    Display = loaddisplay(fname_prefix, Events, iTr, 'StartOn',[0 len], format, experiment.hardware.acquisition.samplingrate);
                    time = (1:len) + Events.StartOn(iTr);
                    time_disp = (1:length(Display))./30 + Events.StartOn(iTr);
                    plot(time, Hand');
                    hold on
                    plot(time_disp, Display/5000)
                    x = repmat(disp_evs(:,2),1,2);
                    plot(x', repmat([-40 20], length(x),1)', 'k--')
                    
                    %plot all the event codes
                    colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1];
                    
                    for i=1:length(goTimes)
                        plot(goTimes(i), 30, 'o', 'color', colors(goID(i)+1,:))
                    end
                    for i=1:length(onTimes)
                        plot(onTimes(i), 40, 's', 'color', colors(onID(i)+1,:))
                    end
                    for i=1:length(acqTimes)
                        plot(acqTimes(i), 50, '*', 'color', colors(acqID(i)+1,:))
                        %plot(acqOffTimes(i), 55, '*', 'color', colors(acqID(i)+1,:))
                    end
                    for i=1:length(touchTimes)
                        plot(touchTimes(i), 60, '>', 'color', colors(1,:))
                    end
                    
                    if ~isnan(seq.rewardOn)
                        plot(seq.rewardOn, 55, '<', 'color', colors(end,:))
                    else
                        plot(seq.errorOn, 55, '<', 'color', [0 0 0])
                    end
                    title(['Trial' num2str(iTr)])
                    %acqPos %#ok<NOPRT>
                    pause
                end
                
            else
                nanseq.acq.time      = nan;
                nanseq.acq.offTime   = nan;
                nanseq.acq.pos       = nan;
                nanseq.acq.targPos   = nan;
                nanseq.acq.targid    = nan;
                nanseq.go.time       = nan;
                nanseq.go.targid     = nan;
                nanseq.on.time       = nan;
                nanseq.on.targid     = nan;
                nanseq.touch.time    = nan;
                nanseq.touch.offTime = nan;
                nanseq.touch.pos     = nan;
                nanseq.touch.targEst = nan;
                nanseq.touch.segid   = nan;
                
                %log reward on, error
                nanseq.rewardOn = nan;
                nanseq.errorOn  = nan;
                
                nanseq.length = nan;
                nanseq.rand   = nan;
                nanseq.index  = nan;
                
                Events.seqEventsParsed(iTr) = nanseq;
            end %if sequence trial
            
        end %trials
        
        %save events file again
        save([fname_prefix '.Events.mat'], 'Events');
        
    end %check for events file
    
end %recordings


%save Trials for the day
saveTrials(day);
