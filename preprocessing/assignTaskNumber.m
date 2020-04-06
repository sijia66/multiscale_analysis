function AdotB_result = assignTaskNumber(Pos_Seq_1_unique,Pos_Seq_1)
% This function calculates the inner product between each target
%  coordnate and defined targets.
%Inputs
%   Pos_Seq_1_unique: number of types of targets * 2 (x,y coordinates)
%   Pos_Seq_1: number of trials * 2 (x,y coordinates)
%Output
%   AdotB_result: number of trials * 1 (assign task)

%vecNorm calculates L2 vector length
A = Pos_Seq_1';
B = Pos_Seq_1_unique';

%vecNorm calculates L2 vector length
B_norm = B ./ repmat(vecnorm(B),size(B,1),1);
A_norm = A ./ repmat(vecnorm(A),size(A,1),1);

AdotB = A_norm' * B_norm;

AdotB(AdotB <= 0.99 ) = 0; % allow for some numerical error
% %check if all trials assigned 
if round(sum(sum(AdotB))) - size(A,2) ~= 0 %
    error('Some trials are not assigned!')
end

%convert to actual task assignment
AdotB_result =  AdotB(:,1);
for ci = 2:size(B,2)
    AdotB_result = AdotB_result + ci * AdotB(:,ci);
end

AdotB_result = uint8(AdotB_result);
end

