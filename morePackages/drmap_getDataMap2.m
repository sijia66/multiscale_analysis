function [D_map, X, Y] = drmap_getDataMap2(D, x, y)

D = D(:);

N = size(D,1);
X = unique(x);
Y = unique(y);

nX = length(X); 
nY = length(Y);


D_map = nan(nX, nY);
for i=1:N
    xid = X==x(i);
    yid = Y==y(i);
    D_map(xid, yid) = D(i);
end