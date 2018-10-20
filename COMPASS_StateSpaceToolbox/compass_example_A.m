
%% Load MSIT data
% data: Ii, Ini, Iin, Yk_1 , Yk
load('MSIT_1.mat');

%% Set behavioral model structure and learning parameters
%  create model
Param = compass_create_state_space(1,0,5,0,[1;1],[1 2],[0 1],[],[]);
% set learning parameters
Iter  = 100;
Param = compass_set_learning_param(Param,Iter,0,1,1,1,1,1,1,2,1);

%% Format the data
% In
In = [Ini Iin  ones(length(Yk),1) Ii Yk_1];
% all data points are valid
valid = ones(length(Yk),1);

%% Run learning with Gamma model
[XSmt,SSmt,Param,XPos,SPos,ML,Yp]=compass_em([2 0],[],In,[],Yk,[],Param,valid);

%% Goodness of fit analysis
[cov_state,cov_obs]=compass_param_covariance_info([2 0],[],In,[],Yk,[],Param,valid,XSmt,SSmt);

%% Extra script
figure(1)
plot(Yk,'LineWidth',2);hold on;
plot(Param.S+Yp,'LineWidth',2);
plot(find(Ii==1),Yk(find(Ii==1)),'o'); 
hold off;
axis tight
legend('RT','Yp','Interfernece Trial')
xlabel('Trial');
ylabel('RT(sec)');

figure(2)
K  = length(Yk);
xm = zeros(K,1);
xb = zeros(K,1);
for i=1:K
    temp=XSmt{i};xm(i)=temp(1);
    temp=SSmt{i};xb(i)=temp(1,1);
end
compass_plot_bound(1,(1:K),xm,(xm-2*sqrt(xb))',(xm+2*sqrt(xb))');
ylabel('x_k');
xlabel('Trial');
axis tight
grid minor


figure(3)
for i=1:Iter
    ml(i)=ML{i}.Total;
end
plot(ml,'LineWidth',2);
ylabel('ML')
xlabel('Iter');

figure(4)
err = Yk - Yp -Param.S;
plot(err,'LineWidth',2);
ylabel('Residual Error')

figure(5)
acf(err,10)

figure(6)
W     =     [Param.Ck(1)   Param.Ck(2)   Param.Dk(3)    Param.Dk(4)    Param.Dk(5)];
W_up  = 1.96*[cov_obs.SE_C(1) cov_obs.SE_C(2) cov_obs.SE_D(3)  cov_obs.SE_D(4)  cov_obs.SE_D(5)];
W_low = -1.96*[cov_obs.SE_C(1) cov_obs.SE_C(2) cov_obs.SE_D(3)  cov_obs.SE_D(4)  cov_obs.SE_D(5)];
errorbar(1:5,W,W_up,W_low,'LineWidth',2);
%xticks([1 2 3 4 5])
ax = gca;
set(ax,'XTick',[1:5])
xticks = get(ax,'XTickLabel');
xticks{1} = 'I_n_2_i';
xticks{2} = 'I_i_2_i';
xticks{3} = 'C_0';
xticks{4} = 'I_i';
xticks{5} = 'Y_k_-_1';
set(ax,'XTickLabel',xticks)
axis tight
