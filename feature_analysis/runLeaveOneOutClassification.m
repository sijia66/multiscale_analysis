function [cat_hat] = runLeaveOneOutClassification(data, categories)

%[cat_hat] = runLeaveOneOutClassification(data, categories)
%
%function to run leave-one-out cross-validated classification. 
%
%inputs: data - matrix (#trials x #features) of data to be classified
%        categories - vector (#trials x 1) of data category labels
%outputs: cat_hat - vector (#trials x 1) of estimated category for each
%                   trial

num_trials = size(data,1);

cat_hat = nan(num_trials,1); %initialize

%loop through trials
for ti = 1:num_trials
    data_row = data(ti,:);
    
    data_temp = data;
    sample_temp = data_temp(ti,:);
    data_temp(ti,:) = [];
    
    cat_temp = categories;
    real_label = cat_temp(ti);
    cat_temp(ti) = [];
    
    cat_hat(ti) = classify(sample_temp,data_temp,cat_temp);
    %Mdl = fitcdiscr(data_temp,cat_temp);
    %cat_hat(ti) = predict(Mdl,data_row);
end

