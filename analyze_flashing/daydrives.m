function [names, bases,recs]=daydrives(day, MonkeyDir)
%
%   [names bases]=daydrives(day, MonkeyDir)
%   checks experiment def files in all recs for the input day. Returns
%   list of microdrives (.name, .name_base fields) that have been recorded
%   on the given day.

global MONKEYDIR

if nargin < 2 || isempty(MonkeyDir)
    MonkeyDir = MONKEYDIR;
end

%get all recording directories
recs = dayrecs(day, MonkeyDir);

%load each exp_def file to get drives
name_list = cell(size(recs));
base_list = cell(size(recs));
for r=1:length(recs)
    
    try
        load([MonkeyDir '/' day '/' recs{r} '/rec' recs{r} '.experiment.mat'], 'experiment')
        name_list{r} = {experiment.hardware.microdrive(:).name};
        
        if isfield(experiment.hardware.microdrive, 'name_base')
            base_list{r} = {experiment.hardware.microdrive(:).name_base};
        end
    catch
        warning(['No experiment definition file found for rec' recs{r}])
    end
end

name_list = cat(2, name_list{:});
base_list = cat(2, base_list{:});

[names, ind] = unique(name_list);

if isempty(base_list)
    bases = [];
else
    bases = base_list(ind);
end
