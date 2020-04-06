function [B,Binfo] = calSlideWindow(A,m,n,A_info)

%this function calculates the sliding window averages for a particular
%matrix, the boundary values are cut off.
% required packages:
%   colfilt in image processing toolbox

%Inputs
%   A:matrix to be averaged
%   m:window size in the first dimension
%   n:window size in the second dimension
%   A_info: values in the x and y directions

%Outputs
%   B: window slided matrix
%   B_info: window averaged axis values

% example 
% A = [1 2 3 4;...
%      5 6 7 9];
% m = 1; 
% n = 2;
% 
% return 
% B = 
% 1.5000    2.5000    3.5000
% 5.5000    6.5000    7.5000
%  B_info = 
%       rowIndex: [1.5000 2.5000 3.5000]
%            colIndex: [1 2]

%Bugs
% colfilt does not seem to support the degenerate case m = n = 1
 
[nRows, nCols ] = size(A);

% rowIndex = 1:nCols;
% colIndex  = 1:nRows;

nRowsStart = max(floor( m / 2 ),1);
nColStart = floor(n / 2);


nRowsOut = nRows - floor(m / 2);
nColsOut = nCols - floor(n / 2);

B = colfilt(A,[m n],'sliding',@mean);
B = B(nRowsStart:nRowsOut,nColStart:nColsOut);


if n == 1
    Binfo.rowIndex = A_info.rowIndex;
else
    rowIndex = colfilt(A_info.rowIndex,[1 n],'sliding',@mean);
    Binfo.rowIndex  = rowIndex(nColStart:nColsOut);
end

if m == 1
    Binfo.colIndex  = A_info.colIndex;
else
    colIndex = colfilt(A_info.colIndex,[1 m],'sliding',@mean);
    Binfo.colIndex  = colIndex(nRowsStart:nRowsOut);
end


end