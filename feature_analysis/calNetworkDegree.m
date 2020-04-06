function  connectDegree = calNetworkDegree(cMatrix,typeOfMatrix)
%this function calculates network degree for each electrode
%Input: 
%       cMatrix: an m by n matrix of binarized network matrix
%       typeOfMatrix: a string
%                     'sym' : for intradrive connections same number of
%                     rows and columsn
%                     'row': calculate connectivity degree along the rows
%                     'col': along the columns
%Output: an m by 1 vector of connectons for each electrode
%        m by 1 vector 
%        n by 1 vector

if strcmp('sym',typeOfMatrix)
    connectDegree = sum(cMatrix,2);
elseif strcmp('row',typeOfMatrix)
    connectDegree = sum(cMatrix,2);
elseif strcmp('row','col')
    connectDegree = sum(cMatrix,1);
else
    disp('Does not recognize the type of the connectivity matrix')
    
end

end

