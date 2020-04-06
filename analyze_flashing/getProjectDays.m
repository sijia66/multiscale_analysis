function subFolders = getProjectDays(dir_name)

files = dir(dir_name);
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

garbage_element = [];

for fi = 1:length(subFolders)
    if contains(subFolders(fi).name,'.')
        garbage_element = [garbage_element fi];
    end
end

subFolders(garbage_element) = [];

subFolders = {subFolders.name};

end