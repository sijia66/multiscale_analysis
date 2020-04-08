% this testing script requires installation of tensorlab
% https://www.tensorlab.net/doc/data.html

%similar to vectorize a matrix
%tensor can be matricized
%the following example demonstrates the example

%first plane
T(:,:,1)= [1 2; 3 4]; 
T(:,:,2)= [5 6;7 8];
T(:,:,3)= [9 10; 11 12];
T(:,:,4)= [13 14;15 16];
size(T)
M_1 = tens2mat(T,[1 2], 3)
size(M_1)

% what happens is that we vectorize dimensions in the first position
% and stack along the columns specified by the second dim
% 
M_2 = tens2mat(T,[1 2 3])

%try different order

M_3 = tens2mat(T, [1 3],2)
M_4 = tens2mat(T, [2 3],1)
M_5 = tens2mat(T, [3 2],1)

%% tensor multiplication

%mode 2 tensor example 
%the matrix is premultiplied by the ccrresponding tensor slice. 
%use the same tensor from above
%mode 1 is simply left multiply each T(:,:,i) slice
U = [1 2 3;4  5 6]';
T_prod_mode1 = tmprod(T,U,1)

%mode 2, the U vector needs to rotate to along the row dimension
%we also need to fold the tensor structure along midpoint of the second
%dimension
T_prod_mode2 = tmprod(T,U,2)

%mode three, if we cut along the third dimension. we have 4 slices into
%T(:,:,i) for i = 1:4
%this operation is simply multiply each row of U_mode3 and add them up
%the first new slice would be i.e. sum of U_mode(1,) * T(:,:,i)
%the same intuition carries
U_mode3 = [1 2 3 4;5 6 7 8];
T_prod_mode3 = tmprod(T,U_mode3,3)

%to illustrate, T_mode
[num_rows, num_cols , num_old_slices] = size(T)
[num_new_slices, num_old_slices] = size(U_mode3)
T_prod_mode3_illus = [];
for j = 1:num_new_slices
    new_slice = zeros(num_rows, num_cols);
    for i = 1:num_old_slices
        old_slice = squeeze(T(:,:,i));
        new_slice = new_slice + old_slice * U_mode3(j,i);
        
    end
    T_prod_mode3_illus(:,:,j) = new_slice;
end

%if all elements are equal, then the logical sum  should equal to their num
%of elements
if sum(T_prod_mode3_illus == T_prod_mode3) == numel(T_prod_mode3)
    disp('the for loop illutration worked')
else
    disp('the for loop illutration failed')
end

%% Kronecker product 
% form a block matrix by compositing each element multiplied the second
% matrix
A = [1 2;3 4];
B = [5 6;7 8];

kron_prod = kron(A,B)

%related column wise Khatri–Rao product 
kr_product = kr(A,B)


%%

X =  T;
figure(1); voxel3(X);
figure(2); surf3(X);
figure(3); slice3(X);


%% experiment with tensor decomp

%generate random tensors 
% "By default, cpd_rnd generates each factor U{n} as randn(size_tens(n),R)" 
size_tens = [7 8 9]; R = 4; 
U = cpd_rnd(size_tens,R);

%generate the associated full tensor
T = cpdgen(U); 

%which is equavalent to
%for each column of U{2} and U{3} 
%we arrange the col of U{3} multiplied by its weight (corresponding element
%in U{2}) as a new column.
%this col is taken the outer product by the corresponding U{1} column
%this repeats for all the columns in U{1}

% M = U{1}*kr(U(end:-1:2)).'; 
% size_tens = cellfun('size',U,1); 
% T = mat2tens(M,size_tens,1);

%% cpd
% now do decomp
[Uhat, output] = cpd(T,R);

semilogy(0:output.Algorithm.iterations,sqrt(2*output.Algorithm.fval)); 
xlabel('iteration'); ylabel('frob(cpdres(T,U))'); grid on;

%these values are not zero possibly PD is not unique
frob(U{1} - Uhat{1})
frob(U{2} - Uhat{2})
frob(U{3} - Uhat{3})

%Uhat needs to be permuted and scaled to match that of the original U
%columns
[relerr,P,D,Uhatps] = cpderr(U,Uhat); 

%can check how many ranks required to capture the origin tensor strucutre
%% try rankest
% this alogoithm calculates  CPD for various ranks and plots the results
rankest(T);



