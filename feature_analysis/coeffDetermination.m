function R_squared = coeffDetermination(X, X_hat, dim)

%R_squared = coeffDetermination(X, X_hat, dim)
%
%computes the coefficient of determination (R^2) between data (X) and its
%estimate (X_hat). Can opperate on a matrix to compute the R^2 for multiple
%components. Calculates variances along the specified dimension 'dim'. 
%
%inputs: X - data matrix, [N M]
%        X_hat - data estimate matrix, [N M]
%        dim   - dimension of matrix to calculate variance along (rows or
%                columns). Default: 1
%outputs: Rsquared - coefficient of determination of X_hat. 
%                    [1 x M] vector if dim = 1
%                    [N x 1] vector if dim = 2

if ~exist('dim', 'var')
    dim = 1;
end

if dim>2
    error('coeffDetermination only defined for vectors or 2-d matrices')
end

%if working along second dimension, flip matrices 
if dim==2
    X = X';
    X_hat = X_hat';
end

ss_res = sum( (X - X_hat).^2,1);

ss_tot = sum( (X - repmat( mean(X,1), size(X,1),1)).^2,1);

R_squared = 1 - ss_res./ss_tot;


%flip back to appropriate orientation
if dim==2
    R_squared = R_squared';
end





