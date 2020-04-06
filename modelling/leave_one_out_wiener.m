function [R_squared_prediction,A_mean] = leave_one_out_wiener(response, data)
%this function performs use ith trial for prediction and the rest trials 
% for training 
% data: number of trials by electrodes by number of times 
%

z = size(response);
z_data = size(data);

trial_seq = 1:z(1);
R_squared_prediction = nan(z(1),z(2));
A_cum = zeros(z_data(2), z(2));

for ti = 1:z(1)

    Y = compress_dim(response(trial_seq~=ti,:,:),2);
    X = compress_dim(data(trial_seq~=ti,:,:),2);
    A = X\Y; %train the model
    
    Y = compress_dim(response(trial_seq==ti,:,:),2);
    X = compress_dim(data(trial_seq==ti,:,:),2);
    
    Y_hat = X*A;
    
    R_squared_prediction(ti,:) = coeffDetermination(Y, Y_hat, 1);
    A_cum = A_cum + A;


end

A_mean = A_cum / z(1);
end