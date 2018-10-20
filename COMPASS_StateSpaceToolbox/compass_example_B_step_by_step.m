%% Load Learning Data
% data: Ii, Ini, Iin, Yk_1 , Yk
load('LEARNING_1.mat');

% Yn - log of reaction time
ind = find(~isnan(Yn));
Yn  = log(Yn(ind)/1000);
Yb  = Yb(ind);
N   = length(Yn);
% Input, 1 xi
In = zeros(N,2);
In(:,1)=1;
In(:,2)=1;
% Input, Ib equal to In
Ib = In;
% Uk, which is zero
Uk = zeros(N,1);
% Valid, which is valid for observed point
Valid = zeros(N,1);
Valid(find(isfinite(Yn)))=1;

%% Set Behavioral Model and Learning Procedure
%  create model
Param = compass_create_state_space(1,1,2,2,eye(1,1),1,1,1,0);
% set learning parameters
Iter  = 250;
Param = compass_set_learning_param(Param,Iter,0,1,0,1,1,1,1,1,0);

Param = compass_set_censor_threshold_proc_mode(Param,log(2),1,1);
 
%% Format the Data
%% Run learning with Gamma model
[XSmt,SSmt,Param,XPos,SPos,ML,YP,YB]=compass_em([1 1],Uk,In,Ib,Yn,Yb,Param,Valid);



%% Extra Script
figure(1)
plot(Yn/1000,'*','LineWidth',2);hold on;
plot(find(Yb),Yb(find(Yb)),'go','LineWidth',2); 
plot(find(Yb==0),Yb(find(Yb==0)),'ro','LineWidth',2); 
grid minor
hold off;
axis tight
legend('RT','True','False','Missing')
xlabel('Trial');
ylabel('RT(sec)/Decision');

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



