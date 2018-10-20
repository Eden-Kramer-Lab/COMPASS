%% Load behavioral data and prepare it for the toolbox
% Data: Yb , Yn
load('model.mat');
% Yn - logarithm of reaction time
Yn = log(Yn(1:end));
Yb = Yb(1:end);
N   = length(Yn);
% Input - 1 xi
In = zeros(N,3);In(:,1)= 1;In(1:2:end,2)= 1;In(:,3)= 1;
% Input, Ib is equal to In
Ib = In;
% Uk, which is all zero
Uk = zeros(N,1);
% valid, which is valid for the observed data point
Valid = ones(N,1); 

%% Build behavioral model and learning procedure
%  create model
Param = compass_create_state_space(2,1,3,3,eye(2,2),[1 2],[0 0],[1 2],[0 0]);
Param.Ek=-1;
% set learning parameters
Iter  = 100;
Param = compass_set_learning_param(Param,Iter,0,1,1,1,1,1,1,1,0);
% define censored point threshold
Param = compass_set_censor_threshold_proc_mode(Param,log(2),1,1);
 
%% Run learning with a mixture of normal & binary
[XSmt,SSmt,Param,XPos,SPos,ML,YP,YB]=compass_em([1 0],Uk,In,Ib,Yn,Yb,Param,Valid);

%% Deviance analysis
%[DEV_C,DEV_D]= compass_deviance([1 1],In,Ib,Yn,Yb,Param,Valid,XSmt,SSmt);

%% Extra Script
figure(1)

K  = length(Yn);
xm = zeros(K,1);
xb = zeros(K,1);
for i=1:K
    temp=XSmt{i};xm(i)=-temp(2)-1;
    temp=SSmt{i};xb(i)=temp(2,2);
end
compass_plot_bound(1,(1:K)*10/K,xm,(xm-sqrt(xb))',(xm+sqrt(xb))');
ylabel('Vigilance');
xlabel('Time (min)');
axis tight
box off
set(gca,'fontsize',28);
set(gca,'FontWeight','Bold');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 27 12])
print('-dpng','f2.png');


figure(3)
for i=1:Iter
    ml(i)=ML{i}.Total;
end
plot(ml,'LineWidth',2);
ylabel('ML')
xlabel('Iter');

% figure(4)
% err = log(Yn/1000) - Yp;
% plot(err,'LineWidth',2);
% ylabel('Residual Err.')
% xlabel('Trial');

% figure(5)
% acf(err,10)

figure(5)
plot(YB,'LineWidth',2)
xlabel('Trial')
ylabel('Prob')
axis tight
grid minor

figure(6)
plot(YP,'LineWidth',2)
xlabel('Trial')
ylabel('RT')
axis tight
grid minor


