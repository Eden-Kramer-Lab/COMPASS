%% Load Model Parameters
load('em_set100');
load('compass_example_var');
%% Initialize Variables
copy_obs_valid = ones(size(obs_valid));
censor_vec = [];
em_list={};
%% Create Full-Likelihood and MAR models
for i = length(th_s):-1:1
    index = find(Yn>th_s(i));
    %Censored Data Treated as Full-Likelihood Model
    obs_valid=copy_obs_valid;
    obs_valid(index)=2;
    Param.censor_time = th_s(i);
    save(['em_set1' num2str(i)],'DISTR','Ib','In','Yn','Yb','Param','obs_valid','Uk')
    file = ['em_set1' num2str(i)];
    em_list = [em_list file];
    censor_vec = [censor_vec length(find(obs_valid==2))];
    %Censored Data Trated as MAR
    obs_valid=copy_obs_valid;
    obs_valid(index)=0;
    save(['em_set2' num2str(i)],'DISTR','Ib','In','Yn','Yb','Param','obs_valid','Uk')
    file = ['em_set2' num2str(i)];
    em_list = [em_list file];
    
end
compass_run_models(em_list) 
%% Initialize Statistic Variables
rmse_a= [];
rmse_b= [];
sx_a  = [];
sx_b  = [];
sy_a  = [];
sy_b  = [];
%% Create Statistics
for i = 1:2:length(th_s)
    %Full-Likelihood
    load(['model_result' num2str(i)])
    sx_a=[sx_a;Param.Wk];
    sy_a=[sy_a;Param.Vk];
    rmse_a=[rmse_a;sqrt(sum((cell2mat(rXSmt)-x').^2))];
    %MAR
    load(['model_result' num2str(i+1)])
    sx_b=[sx_b;Param.Wk];
    sy_b=[sy_b;Param.Vk];
    rmse_b=[rmse_b;sqrt(sum((cell2mat(rXSmt)-x').^2))];
   
end
%% Save Statistics
save('compass_example_E_variables','rmse_a','rmse_b','sx_a','sx_b','sy_a','sy_b','censor_vec')
%% Create Graphs of Observed vs. Estimated Statistics
figure
subplot(4,1,1)
plot(rmse_a,'r','LineWidth',2);hold on;plot(rmse_b,'b','LineWidth',2);title('RMSE Error');hold off;
subplot(4,1,2)
plot(sx_a,'r','LineWidth',2);hold on;plot(sx_b,'b','LineWidth',2);title('State Noise Paramater Estimate');hold on;
plot(sx*ones(size(sx_a)));hold off;
subplot(4,1,3)
plot(sy_a,'r','LineWidth',2);hold on;plot(sy_b,'b','LineWidth',2);title('Observation Noise Paramater Estimate');hold on;
plot(sy*ones(size(sy_a)));hold off;
legend('Likelihood','MAR','True Parameter');
subplot(4,1,4)
plot(censor_vec,'b','LineWidth',2,'Color','r'); title('Number of Censored Points')
xlim([1 20])
xlabel('Model No.')

figure
plot(x)
hold on
plot(y)
legend('Estimate','True Yn Value')
xlabel('Iteration No.')
ylabel('Time (s)')
title('True Yn Value over Estimate')
