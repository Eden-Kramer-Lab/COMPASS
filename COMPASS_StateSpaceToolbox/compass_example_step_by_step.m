%% Load behavioral data and prepare it for the toolbox
% Data: Yb (decision) , Yn (reaction time)
load('LEARNING_1.mat');
% Yn - logarithm of reaction time
Yn  = log(Yn/1000);
% Data length
N   = length(Yn);
% Input - 1 xi
In = zeros(N,2);In(:,1)= 1; In(:,2)= 1;
% Input, Ib is equal to In
Ib = In;
% Uk, which is all zero
Uk = zeros(N,1);
% valid, which is valid for the observed data point
ind = find(isnan(Yn)); 
Valid = ones(N,1);
Valid(ind)=0;

%% Data Visualization
figure(1)
plot(Yn,'LineWidth',2);hold on;plot(ind,-2,'kx','LineWidth',2);
ylabel('Reaction Time')
yyaxis right
plot(Yb,'o','LineWidth',2);
ylabel('Decision Choice (0/1)')
xlabel('Trial Index')
axis tight
title('Behavioral Signal')

%% Build behavioral model and learning procedure
%  create model
Param = compass_create_state_space(1,1,2,2,1,1,1,1,0);
% set learning parameters
Iter  = 250;
Param = compass_set_learning_param(Param,Iter,0,1,0,1,1,1,1,1,0);
% define censored point threshold
Param = compass_set_censor_threshold_proc_mode(Param,log(2),1,1);
 
%% Run learning with a mixture of normal & binary
% note, we ran compass_em once with default setting for parameters and then
% Initialize parameters using GLM method. 
Param.Wk = 0.2;
Param.Ck = 0.1;
Param.Dk = [0 -1.66];
Param.Vk = 0.05;
Parm.Fk = [0 -.15];
[XSmt,SSmt,Param,XPos,SPos,ML,YP,YB]=compass_em([1 1],Uk,In,Ib,Yn,Yb,Param,Valid);

figure(2)
K  = length(Yn);
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
title('Learning State')


figure(3)
for i=1:Iter
    ml(i)=ML{i}.Total;
end
%subplot(141)
plot(ml,'LineWidth',3);
axis tight
%ylabel('Total ML')
ylabel('ML')
xlabel('Iter');
for i=1:Iter
    ml(i)=ML{i}.State;
end
title('Maximum Likelihood Curve')
% subplot(142)
% plot(ml,'LineWidth',3);
% axis tight
% ylabel('State ML')
% xlabel('Iter');
% for i=1:Iter
%     ml(i)=ML{i}.ObsrvNormGamma;
% end
% subplot(143)
% plot(ml,'LineWidth',3);
% axis tight
% ylabel('Cont. ML')
% xlabel('Iter');
% for i=1:Iter
%     ml(i)=ML{i}.ObsrvBern;
% end
% subplot(144)
% plot(ml,'LineWidth',3);
% axis tight
% ylabel('Disc. ML')
% xlabel('Iter');

figure(4)
subplot(121)
plot(YP,'LineWidth',2);hold on;plot(Yn,'LineWidth',2) ;hold off;
xlabel('Trial')
ylabel('Reaction Time')
axis tight
title('RT plus its prediction')

subplot(122)
Yb(ind)=nan;
plot(YB,'LineWidth',2);hold on;plot(Yb,'o','LineWidth',2) ;hold off;
xlabel('Trial')
ylabel('Decision Choice')
axis tight
title('Decision plus Probability of Correct Choice')


%% Deviance analysis
[DEV_C,DEV_D]= compass_deviance([1 1],In,Ib,Yn,Yb,Param,Valid,XSmt,SSmt);

% Param.Wk = 0.2;
% Param.Ck = 0.1;
% Param.Dk = [0 -1.66];
% Param.Vk = 0.05;
% Param.Ek = 0;
% Parm.Fk  = [0 -.15];
% [XSmt,SSmt,Param,XPos,SPos,ML,YP,YB]=compass_em([1 1],Uk,In,Ib,Yn,Yb,Param,Valid);
% [DEV_C1,DEV_D1]= compass_deviance([1 1],In,Ib,Yn,Yb,Param,Valid,XSmt,SSmt);
% 
% Param.Wk = 0.2;
% Param.Ck = 0.1;
% Param.Dk = [0 -1.66];
% Param.Vk = 0.05;
% Param.Ek = 1;
% Parm.Fk  = [0 -.15];
% [XSmt,SSmt,Param,XPos,SPos,ML,YP,YB]=compass_em([0 1],Uk,In,Ib,Yn,Yb,Param,Valid);
% [DEV_C2,DEV_D2]= compass_deviance([1 1],In,Ib,Yn,Yb,Param,Valid,XSmt,SSmt);


%% Covariance analysis
[COV_X,COV_C,COV_D]=compass_param_covariance_info([1 1],Uk,In,Ib,Yn,Yb,Param,Valid,XSmt,SSmt);

Ps = [Param.Wk          Param.Ck     Param.Dk(2)     Param.Vk     Param.Fk(2)];
Ss = [COV_X{1}.SE_W     COV_C.SE_C   COV_C.SE_D(2)   COV_C.SE_V   COV_D.SE_F(2)];

figure(5)
errorbar(1:5,Ps,2*Ss,'.','LineWidth',3)
grid on
axis tight
set(gca,'XTick',[1 2 3 4 5])
set(gca,'XTickLabel',{'Sv','b1','b0','Sw','c0'})
%xticks([1 2 3 4 5])
xticklabels({'Sw','b1','b0','Sv','c0'})
