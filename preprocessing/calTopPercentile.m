function perCols = calTopPercentile(cols, perCentile)
%this function outputs a column vector containing top percentile elements
%in column cols
%Inputs:
%   cols: input column vector n by 1
%   perCentile: the specified percentile
%Output
%   perCols: a x by 1 column vector by 1,  x is the top % elements


if ~iscolumn(cols)
    disp('the input is not a column vector, return')
    perCols = [];
    return
end

prctileValue = prctile(cols,perCentile);
%select the values
perCols = find(cols >= prctileValue);

end