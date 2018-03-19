close all
clear all
%% Load Reinforcement Learning Data
% data: In, RESP, dGain, Type, List, Gain
% There are 160 trials in total
% In:    160x4, each column corresponds to one of four different states
% Resp:  160x1, it is the action- there are two actions (0 or 1)
% dGain: 160x1, it is the reward received per each trial
% Type:  160x1, it indicates loss-win trials (1 is for the win trials, and -1 is for the loss trials)
% List:  160x1, it show the trial state - it matches to In
% Gain:  160x1, is the comulative gain
%load('example_b.mat');
load('extra_example_b.mat')
In = zeros(160,4);
In(find(List==1),1)=1;
In(find(List==2),2)=1;
In(find(List==3),3)=1;
In(find(List==4),4)=1;


%% Set Behavioral Model
% create model
% - 2 state variables per state plus one for action
% - six input to the state transition process (check the paper for details)
% - 14 input in the obsevation process (check the paper for details)
%   * input 1 and 2 are for state 1 - input 1 is all zeros, and input 2 is  all one
%   * input 3 and 4 are for state 2 - input 3 is all zeros, and input 4 is  all one    
%   * input 5 and 6 are for state 3 - input 5 is all zeros, and input 6 is  all one    
%   * input 7 and 8 are for state 4 - input 7 is all zeros, and input 8 is  all one     
% - input 9 is always one linked to Q state variable function
% - input 10 is always one linked to baseline of Q function 
% - input 11 is one on state 1 trials 
% - input 12 is one on state 2 trials 
% - input 13 is one on state 3 trials 
% - input 14 is one on state 4 trials 
Param = compass_create_state_space(9,4,0,9,eye(9,9),[],[],[1 2 3 4 5 6 7 8 9],[0 0 0 0 0 0 0 0 0]);
RATE_A = 0.1;
RATE_B = 0.1;
RATE_C = 0.1;
% set initial value for Wv,s1 to RATE_A
Param.Ak(1,1)=  1-RATE_A; 
Param.Ak(2,1)=  1;          Param.Ak(2,2)= -RATE_B;
% set initial value for Wv,s2 to RATE_A
Param.Ak(3,3)=  1-RATE_A;
Param.Ak(4,3)=  1;          Param.Ak(4,4) = -RATE_B;
% set initial value for Wv,s3 to RATE_A
Param.Ak(5,5)= 1-RATE_A;
Param.Ak(6,5)= 1;           Param.Ak(6,6)= - RATE_B;
% set initial value for Wv,s4 to RATE_A
Param.Ak(7,7)= 1-RATE_A;
Param.Ak(8,7)= 1;           Param.Ak(8,8)=  - RATE_B;
% set initial value for Wq,a1 to 1
Param.Ak(9,9)= 1-RATE_C;
% set initial value of B
Param.Bk(1,1)   = RATE_A;
Param.Bk(2,1)   = RATE_B;
Param.Bk(3,2)   = RATE_A;
Param.Bk(4,2)   = RATE_B;
Param.Bk(5,3)   = RATE_A;
Param.Bk(6,3)   = RATE_B;
Param.Bk(7,4)   = RATE_A;
Param.Bk(8,4)   = RATE_A;
Param.Bk(9,1:4) = RATE_C;
% set Param.Wk
Param.Wk = 0.01 * eye(9,9);
% Set Param.Ek
Param.Ek(1) = 0;Param.Ek(3) = 0;Param.Ek(5) = 0;Param.Ek(7) = 0;


%% Set Learning Procedure
Iter  = 1000;
Param = compass_set_learning_param(Param,Iter,3,1,1,1,1,1,0,2,0);

%% Format the Data
% Yb
Yb = Resp;
% Ib
Ib = zeros(length(Yb),9);
Ib(find(In(:,1)),2) = 1;
Ib(find(In(:,2)),4) = 1;
Ib(find(In(:,3)),6) = 1;
Ib(find(In(:,4)),8) = 1;
Ib(:,9) = 1; 
% all data points are valid
valid = ones(length(Yb),1);

% Uk - is the reward
Uk = zeros(length(Yb),4);
Uk(:,1)= dGain.*Ib(:,2);
Uk(:,2)= dGain.*Ib(:,4);
Uk(:,3)= dGain.*Ib(:,6);
Uk(:,4)= dGain.*Ib(:,8);



%% Run learning with Gamma model
[XSmt,SSmt,Param,XPos,SPos,ML,~,Pb]=compass_em([0 1],Uk,[],Ib,[],Yb,Param,valid);


%% Plot Modeling Result
figure(1)
subplot(1,2,1)
imagesc(Param.Ak);
colorbar
title('Ak');
title('Ak Coefficients')

subplot(1,2,2)
imagesc(Param.Bk);
colorbar
title('Bk');
title('Bk Coefficients')

figure(2)
ml = [];
for i=1:Iter
    ml(i)=ML{i}.Total;
end
plot(ml,'LineWidth',2);
ylabel('ML')
xlabel('Iter');
title('Maximum Likelihood Curve')

figure(3)
plot(Yb,'b--','LineWidth',0.1);hold on;
plot(Yb,'bo','LineWidth',1);
plot(Pb,'r','LineWidth',2);hold off;
ylabel('Yb/Pb')
xlabel('Trial');
title('Decision and Probability of Correct Choice')


figure(4)
K  = length(Yb);
for s=1:9
    if s<9
        d = floor(s/2);
        if mod(s,2)
            subplot(2,5,1+d)
        else
            subplot(2,5,5+d)
        end
    else
             subplot(2,5,10)
    end
    xm = zeros(K,1);
    xb = zeros(K,1);
    for i=1:K
        temp= XSmt{i}; xm(i)=temp(s);
        temp= SSmt{i}; xb(i)=temp(s,s);
    end
    compass_plot_bound(1,(1:K),xm,(xm-2*sqrt(xb))',(xm+2*sqrt(xb))');
    ylabel(['x_'  num2str(s)]);
    xlabel('Trial');
    axis tight
    grid minor
end

% Here, We show Xa*w-Xs
figure(5)
for s=1:4
    Xa=[];
    Xs=[];
    for i=1:K
        temp= XSmt{i}; 
        xm(i)=temp(9)*Param.Ek(end)-temp(s*2);
        temp= SSmt{i}; 
        xb(i)=temp(9,9)*Param.Ek(end)*Param.Ek(end)+temp(2*s,2*s);
    end
    subplot(1,4,s)
    compass_plot_bound(1,(1:K),xm,(xm-2*sqrt(xb))',(xm+2*sqrt(xb))');
    hold on;
    plot(1:K,zeros(K,1),'r')
    grid on
    grid minor
    axis tight
    xlabel('trial')
    ylabel(['Xa *  w - Xs_'  num2str(s)])
    ylim([-3 2.5])
end


figure(6)
plot(dGain,'LineWidth',2);
xlabel('Trial')
ylabel('Gain')
grid minor
title('Gain Over Trials')

figure(7)
ind=find(List==1);plot(ind,-0.05+Resp(ind),'b*','LineWidth',2,'MarkerSize',12);hold on;
ind=find(List==2);plot(ind,-0.1+Resp(ind),'go','LineWidth',2,'MarkerSize',12);hold on;
ind=find(List==3);plot(ind,0.05+Resp(ind),'r+','LineWidth',2,'MarkerSize',12);hold on;
ind=find(List==4);plot(ind,0.1+Resp(ind),'mx','LineWidth',2,'MarkerSize',12);hold off;
legend('Reward Item 1','Reward Item 2','Loss Item 1','Loss Item 2','location','best')
xlabel('Trial');
ylabel('Correct/incorrect choice')





