function [scatterMap,mapInfo] = scatterToMatrix(P_chamber,values)
% this function maps a list of coordinates into a matrix
%Inputs:
%       P_chamber: n by 2 matrix with coordinates
%       vallues: n by 1 column vector for the values
%Output:
%       scatterMap: mapped values     
%       mapInfo: a structure with pitch values and range
%Example Usage:
% [map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,connectDegree);
% imagesc(map)
% axis xy
% %       only label x axis
% [x_num,y_num] = size(map);
% xticks = 1:x_num;
% xlabels = ((0:x_num-1)) * mapInfo.x_pitch + mapInfo.x_zero;
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
% xlabel('mm')
%
%OR
%[map,mapInfo]= scatterToMatrix(P_chamber_ECOG_3,connectDegree);
%imagesc(map)
%scatterToMatrix_plot_label
%
%Author: Si Jia Li Aug 21/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%number of coordinates 
n_coor = size(P_chamber,1);

%do some check up
if n_coor ~= size(values, 1)
    error('The coordinate matrix and the value vector don''t have the same number of rows')
end

% coor_max = max(P_chamber);
% coor_min = min(P_chamber);

x_coors = P_chamber(:,1);
y_coors = P_chamber(:,2);

x_coors_sorted = sort(unique(x_coors));
y_coors_sorted = sort(unique(y_coors));

x_coors_sorted_diff = diff(x_coors_sorted);
y_coors_sorted_diff = diff(y_coors_sorted);


if ~all(round(x_coors_sorted_diff,4) == round(x_coors_sorted_diff(1),4)) %round to circumvent the numerical precision issue
    error('The pitch between electrodes is not the same in the first dimension')
end

if ~all(round(y_coors_sorted_diff,4) == round(y_coors_sorted_diff(1),4))
    error('The pitch between electrodes is not the same in the second dimension')
end

x_pitch = x_coors_sorted_diff(1);
y_pitch = y_coors_sorted_diff(1);

%now we can generate the map
scatterMap = NaN(length(x_coors_sorted),length(x_coors_sorted));

%assign the values 
x_zero = x_coors_sorted(1);
y_zero = y_coors_sorted(1);

for ni = 1:n_coor
    x_index = int8((x_coors(ni) - x_zero) / x_pitch + 1);
    y_index = int8((y_coors(ni) - y_zero) / y_pitch + 1);
    scatterMap(x_index,y_index) = values(ni);
end

%collect information for reporting
mapInfo.x_pitch = x_pitch;
mapInfo.y_pitch = y_pitch;
mapInfo.x_zero = x_zero;
mapInfo.y_zero = y_zero;

end