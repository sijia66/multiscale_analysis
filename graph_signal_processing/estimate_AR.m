%% esitmate AR coefficients

% try matlab example
y = sin((1:300)') + 0.5*randn(300,1);
y = iddata(y);
sys_b = ar(y,4,'burg');



compare(y,sys_b,5)

%get the aic information
value = aic(sys_b)

%% select max model order

trial_1 = 1;
electrode_1 = 2;

Ts = 0.2; %sampling time
model_order = 60; %from the paper 30

data_trial = sq(data_twoDrives{1,2}(trial_1,electrode_1,:));
data_trial_struct = iddata(data_trial,[],Ts);


aic_models = zeros(model_order,1);
acc_models = zeros(model_order,1);

for model_i = 1:model_order
    %fit to the model
    sys_test = ar(y,model_i);
    
    %save AIC
    aic_models(model_i) = aic(sys_test);
    
    %save accuracy
    acc_models(model_i) = sys_test.Report.Fit.FitPercent;
    
    if mod(model_i,10) == 0
        fprintf('Finished: %d models\n',model_i)
    end
end

subplot(2,1,1)
plot(aic_models)
title('Comparing AIC for AR models of Increasing Model Order')
xlabel('AR Model Order')
ylabel('AIC')

subplot(2,1,2)
plot(acc_models)
title('Comparing Accuracy for AR models of Increasing Model Order')
xlabel('AR Model Order')
ylabel('Accuracy (%)')




