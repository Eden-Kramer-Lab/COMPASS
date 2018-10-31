%% Load Parameters
load('compass_example_var');
load('compass_example_E_variables')
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