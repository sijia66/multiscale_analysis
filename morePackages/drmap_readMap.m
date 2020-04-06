function [ids, labels] = drmap_readMap(file, sheet, in, out)

%reads excel file ("file", sheet # "sheet") with mapping info. Collects 'in' column and 'out'
%columns. Returns that data in ids, and labels with names of columns


[num, txt] = xlsread(file, sheet);

labels = cat(2, in, out);

%find columns. Assure order is same as 'labels'
col = nan(1, size(labels,2));
for i=1:size(labels,2);
    
col(i) = find(ismember(lower(txt(1,:)), labels{i}));
end

%grab out relevant data
ids = num(:, col);