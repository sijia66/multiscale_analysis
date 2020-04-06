function mi_x_y = calMI_xy(x,y,xNeglected, yNeglected)
%this function calculates mutual information between two matrices x and y
%   Inputs:
%   x, number of trials * number of electrodes N1,
%   y, number of trials * number of electrodes N2
%   x_neglected: a vector containing the x indices to be neglected:mi
%   y_neglated: a vector containing the y indices to be neglected
%
%   Required packages:
%   MI_toolbox
%
%   output:
%   mi_x_y: N1 by N2
%
%   To be continued:
%   check if electrode has gaussian distribution
%

%separate values into discrete bins
X = ceil(x);
Y = ceil(y);

[xR,xC] = size(x);
[yR,yC] = size(y);

if xR ~= yR
    error('x and y do not have the same number of trials')
end

mi_x_y = zeros(xC,yC);

for ei = 1:xC
    if ismember(ei,xNeglected)
        fprintf('Neglect X: %d\n',ei)
        continue
    end
    for ej = 1:yC
        if ismember(ej,yNeglected) %check if this electrode should be neglected

            mi_x_y(ei,ej) = 0;
            continue
        end
       
        
        mi_x_y(ei,ej) = mi(X(:,ei),Y(:,ej));
    end
end

end

