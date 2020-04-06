function Session = makeDriveDatabase(drive_base, date_range, monkeydir)

% Session = makeDriveDatabase(name_base, date_range, monkeydir)
%
%This function looks through all day directories within date_range (or all in the directory if date-range not specified)
%and makes a 'Drive Database'. The database is a Sessions cell with the day and recordings where the specified drive was recorded
%Note that this searches by 'name_base' and not for a specific drive.
%
%input: drive_base  - string, name_base of drive to find.
%       date_range - cell of strings, {'YYMMDD', 'YYMMDD'} specifying the start and end date ranges to include. 
%                    To include all days, either set to 'all' or {}. If not specified, defaults to 'all'.
%       monkeydir  - string w/ project directory. If not input or empty, uses the global MONKEYDIR
%     
%output: Session - session cell w/ session info for a given drive
%                  Cell format for each session:
%                  {DAY, {recs}, DRIVE_NAME, DRIVE_BASE, IMPLANT_DATE, DRIVE_ID, DRIVE_THETA, DRIVE_X, DRIVE_Y, SESS_#}
%                  A session is created for each unique drive (i.e. drive name is different, but same drive base) and for each day.

global MONKEYDIR

if ~exist('date_range', 'var') || isempty(date_range)
    date_range = 'all';
else
    if ~iscell(date_range) || length(date_range)~=2
        error('date_range must be a cell of length 2 w/ start and end dates (strings, ''YYMMDD'' format)')
    end
end
if ~exist('monkeydir', 'var') || isempty(monkeydir)
    monkeydir = MONKEYDIR;
end


%get all day directories w/in the project
days = getProjectDays(monkeydir);

%trim down to days w/in range
if ~isequal(date_range, 'all')
    
    minDate = str2double(date_range{1});
    maxDate = str2double(date_range{2});
    
    if ~(minDate <= maxDate)
        error('max date must be larger than minimum date')
    end
    
    keepday = false(size(days));
    for iDay = 1:length(days)
        keepday(iDay) = str2double(days{iDay}) >= minDate & str2double(days{iDay}) <= maxDate;
    end 
    
    days = days(keepday);
    
end

%now for each day, find if the drive of interest is recorded
Session = cell(length(days)*2,1);
cnt = 0;
for iDay = 1:length(days)
    
    %find drives recorded that day and recs associated w/ each drive
    [names, bases, recs] = daydrives(days{iDay}, monkeydir); 
    
    
    %if one of the recorded drives matches the name_base requested for this
    %database, add an entry to the database
    if ismember(drive_base, bases)
       
        %indices where drive base matches requested name_base
        %note that multiple might match because daydrives groups output by
        %unique 'name' not unique 'base'. We also want to group/split our
        %sessions by unique name, not unique name_base
        idx = find( ismember(bases, drive_base));
        
        for iDrive = idx(:)
            
        
        %look at recordings where this particular drive was recorded and get all metadata
        [~, ~, day_impDate_list, day_id_list, day_theta_list, day_x_list,day_y_list] =...
            expDriveInfo(days{iDay}, recs, drive_base, monkeydir);
        
        impDate = unique(day_impDate_list); 
        id      = unique(day_id_list);
        theta   = unique(day_theta_list);
        x       = unique(day_x_list);
        y       = unique(day_y_list);
        
        if length(impDate)>1 || length(id)>1 || length(theta)>1 || length(x)>1 || length(y)>1
            error('something is wrong with drive indexing.') 
        end
        
        Session{cnt+1} = {days{iDay}, recs, names{iDrive}, bases{iDrive}, impDate{:}, id{:}, theta, x, y, cnt+1}; %#ok<CCAT>
        
           cnt = cnt+1;
        end %across drives to add
        
    end %if drive to be added to database
    
end %across days

Session = Session(1:cnt);
