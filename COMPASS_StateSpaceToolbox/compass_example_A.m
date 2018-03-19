%% Load behavioral data and prepare it for the toolbox
% Data: Yb , Yn
load('example_a.mat');
% Yn - logarithm of reaction time

ind = find(~isnan(Yn));
Yo  = Yn(ind)/1000;
Yn  = log(Yn(ind)/1000);
% Yb - decision (0/1)
Yb  = Yb(ind); N   = length(Yn);
% Input - 1 xi
In = zeros(N,2);In(:,1)= 1; In(:,2)= 1;
% Input, Ib is equal to In
Ib = In;
% Uk, which is all zero
Uk = zeros(N,1);
% valid, which is valid for the observed data point
Valid = ones(N,1);  

%% Build behavioral model and learning procedure
%  create model
Param = compass_create_state_space(1,1,2,2,eye(1,1),1,1,1,0);
% set learning parameters
Iter  = 250;
Param = compass_set_learning_param(Param,Iter,0,1,0,1,1,1,1,1,0);
% define censored point threshold
Param = compass_set_censor_threshold_proc_mode(Param,log(2),1,1);
 
%% Run learning with a mixture of normal & binary
[XSmt,SSmt,Param,XPos,SPos,ML,YP,YB]=compass_em([1 1],Uk,In,Ib,Yn,Yb,Param,Valid);

%% Deviance analysis
[DEV_C,DEV_D]= compass_deviance([1 1],In,Ib,Yn,Yb,Param,Valid,XSmt,SSmt);

%% Extra Script
figure(1)
plot(find(Yb),Yb(find(Yb)),'go','MarkerSize',8,'LineWidth',3);hold on 
plot(find(Yb==0),Yb(find(Yb==0)),'ro','MarkerSize',8,'LineWidth',3); 
ylabel('Decision')
set(gca,'YTick',[0 1])
set(gca,'YTickLabel',{'Incorrect','Correct'})
set(gca,'YTickLabelRotation',90)
hold on
yyaxis right
plot(1:length(Yn),Yo,'LineWidth',3);hold on;
ylabel('Reaction Time (sec)')
xlabel('Trial Index')
hold off
title('Behavioral Signal')


figure(2)
plot(YB,'LineWidth',3);hold on;
hold on;
load('prerau.mat');
plot(bmode(2:end),'LineWidth',3);hold off
ylabel('Probability of Picking Correct Choice')
xlabel('Trial Index');
hold off
legend('COMPASS','Prerau Code')
title('Learnign State')

figure(3)
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
title('Learnign State plus Its Confidence Interval')


figure(5)
for i=1:Iter
    ml(i)=ML{i}.Total;
end
plot(ml,'LineWidth',2);
ylabel('ML')
xlabel('Iter');
title('Maximum Likelihood Curve')





