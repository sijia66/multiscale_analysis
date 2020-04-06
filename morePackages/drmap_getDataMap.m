function [D_map, D_map_xy, X, Y] = drmap_getDataMap(D, x, y, r, c)

D = D(:);

N = size(D,1);
R = unique(r);
C = unique(c);
if min(R) == 0
    r = r+1;
end
if min(C)==0
    c = c+1;
end
X = unique(x);
Y = unique(y);

D_map = nan(length(R), length(C));
x_map = nan(length(R), length(C));
y_map = nan(length(R), length(C));
for i=1:N    
    D_map(r(i), c(i)) = D(i);
    x_map(r(i), c(i)) = x(i);
    y_map(r(i), c(i)) = y(i);
end

%figure out if x/y correspond to row or columns to make D_map_xy
D_map_xy = D_map;

tmp1 = mode(diff(x_map, 1, 1),1);
tmp2 = mode(diff(y_map, 1, 2),2);
if nansum(tmp1)==0 & nansum(tmp2)==0
    %x = row, y = col. No need to change
    x_dim = 1;
    y_dim = 2;
    
    x_list = mode(x_map, x_dim);
    y_list = mode(y_map, y_dim);
    
    x_dir = x_list(end) - x_list(1);
    y_dir = y_list(end) - y_list(1);
    if x_dir < 0
        D_map_xy = flip(D_map_xy, x_dim);
    end
    if y_dir < 0
        D_map_xy = flip(D_map_xy, y_dim);
    end
    
    %x and y are set appropriately so no need to transpose
else
    %x = col, y = row. swap the map around
    x_dim = 2;
    y_dim = 1;

    x_list = mode(x_map, x_dim);
    y_list = mode(y_map, y_dim);
    
    x_dir = x_list(end) - x_list(1);
    y_dir = y_list(end) - y_list(1);
    if x_dir < 0
        D_map_xy = flip(D_map_xy, x_dim);
    end
    if y_dir < 0
        D_map_xy = flip(D_map_xy, y_dim);
    end
    
    %transpose since x and y are flipped w.r.t. row/col
    D_map_xy = D_map_xy';
end

    