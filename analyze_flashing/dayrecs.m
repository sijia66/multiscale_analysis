function recs=dayrecs(day, MonkeyDir)
%
%   recs=dayrecs(day, MonkeyDir)
%
global MONKEYDIR

if nargin < 2 || isempty(MonkeyDir)
    MonkeyDir = MONKEYDIR;
end

if iscell(day), day = day{1}; end
tmp = dir([MonkeyDir '/' day '/0*']);
for i = 1:8 %only go up to 8 b/c 901+ reserved for 'night' recordings
 tmp = [tmp;dir([MonkeyDir '/' day '/' num2str(i) '*'])];
end

tmp = tmp([tmp(:).isdir]); %only keep directories. 

if ~isempty(tmp)
    [recs{1:length(tmp)}] = deal(tmp.name);
else
    recs = {};
end
