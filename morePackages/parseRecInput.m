function RECS = parseRecInput(recs)

%RECS = parseRecInput(recs)
%This function takes 'recs' input of various forms and puts it into a
%standardized format (cell of strings w/ list of 'NNN' rec ids). This allows each analyze function to take flexible
%'rec' inputs from user but then operate on recs in a standard format.
%
%input: recs - accepts following possible formats
%                  cell of strings ('NNN' or 'N' format)
%                  string ('NNN' or 'N' format)
%                  numeric (vector of numbers). Note: this assumes you're
%                      referencing the REC # (i.e. recs = [1 3] returns RECS =
%                      {'001', '003'} regardless of the list of recs in a given day.
%

RECS = {};

%if it's a cell, loop through elements to check each
if iscell(recs)
    
    
    for r=1:length(recs)
        
        %if string, double-check using 'NNN' format
        if ischar(recs{r})
            if length(recs{r})~=3
                recnum = str2double(recs{r});
                prenum = [];
                if recnum < 10
                    prenum = '00';
                elseif recnum < 100
                    prenum = '0';
                end
                
                RECS{r} = [prenum recs{r}];
            else
                RECS{r} = recs{r};
            end
            
            %if a number, convert to 'NNN' string
        elseif isnumeric(recs{r})
            
            recnum = recs{r};
            prenum = [];
            if recnum < 10
                prenum = '00';
            elseif recnum < 100
                prenum = '0';
            end
            
            RECS{r} = [prenum num2str(recs{r})];
        else
            warning('parseRecInput can only deal with numeric or string rec formats.')
        end
    end
    
    %if not a cell
else
    %number/vector --> string
    if isnumeric(recs)
        for r=1:length(recs)
            recnum = recs(r);
            prenum = [];
            if recnum < 10
                prenum = '00';
            elseif recnum < 100
                prenum = '0';
            end
            
            RECS{r} = [prenum num2str(recs(r))];
        end
        
        %if its a string, double-check using 'NNN' format and pack into a
        %cell.
    elseif ischar(recs)
        
        if length(recs)~=3
            recnum = str2double(recs);
            prenum = [];
            if recnum < 10
                prenum = '00';
            elseif recnum < 100
                prenum = '0';
            end
            
            RECS{1} = [prenum recs];
        else
            RECS{1} = recs;
        end
    else
        warning('parseRecInput can only deal with numeric or string rec formats.')
    end
end
