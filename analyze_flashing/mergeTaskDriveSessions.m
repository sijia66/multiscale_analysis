function Session = mergeTaskDriveSessions(taskSession, driveSession)

%this function takes a 'task'-based Session list and a 'drive-based'
%Session list and merges the information together. The output merged
%Session will be split up into {day, recording} pairs where the drive was
%also recorded. In cases where multiple unique drives were recorded within
%a day (e.g. if an ecog array was removed and re-implanted and recorded
%from pre- and post-replacement), Sessions will divide recordings for that
%day into separate Sessions for each drive. 
%
%input: taskSession - cell w/ Session info for a task. Format of each Session{ind} = {day, recs, ..., ind}
%       driveSession - cell w/ Session infor for drive. 
%                      Format of each Session{ind} = {DAY, {recs}, DRIVE_NAME, DRIVE_BASE, IMPLANT_DATE, DRIVE_ID, DRIVE_THETA, DRIVE_X, DRIVE_Y, ind}
%
%output: Session - cell w/ merged Session info. Format:
%                  Session{ind} = {DAY, {recs}, [any task session info...], DRIVE_NAME, DRIVE_BASE, IMPLANT_DATE, DRIVE_ID, DRIVE_THETA, DRIVE_X, DRIVE_Y, ind}
%


%turn taskSessions and driveSessions into #sessions x N cells -- makes
%getting info out easier w/out loops
taskSession = cat(1, taskSession{:});
driveSession = cat(1, driveSession{:});

[~, N_taskDat]   = size(taskSession);  % # task data in a session
[~, N_driveDat] = size(driveSession); % # drive data in a session



%get day, rec info forTask and Drive sessions
day_task  = taskSession(:,1);
day_drive = driveSession(:,1);

recs_task  = taskSession(:,2);
recs_drive = driveSession(:,2);


%get days shared between the two
commonDays = intersect(day_task, day_drive);
nCommDays = length(commonDays);

%loop to add common days into merged Session
Session = cell(nCommDays*2,1);
cnt = 0;
for iDay = 1:nCommDays
    
    %indices within task and drive sessions
    idx_task  = find(ismember(day_task, commonDays{iDay}));
    idx_drive = find(ismember(day_drive, commonDays{iDay}));
    
    %loop through all task sessions
    for iTask = idx_task(:)
        
        for iDrive = idx_drive(:)
            
            %get common recordings
            commonRecs = intersect(recs_task{iTask}, recs_drive{iDrive});
            
            Session{cnt+1} = {commonDays{iDay}, commonRecs, taskSession{iTask,3:N_taskDat-1}, driveSession{iDrive, 3:N_driveDat-1}, cnt+1};
            
            cnt = cnt+1;
        end
    end
end

Session = Session(1:cnt);
