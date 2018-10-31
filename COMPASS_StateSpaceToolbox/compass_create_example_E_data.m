clear all
%% Set Initial Parameters
sx = 0.01;
sy = 0.055;
K  = 200;
a=0.99;
x0 = sqrt(sx)*randn();
%Create estimate of Yn
for k=1:K
    if k==1
        x(k)=a*x0 + sqrt(sx)*randn();
    else
        x(k)=a*x(k-1) + sqrt(sx)*randn();
    end
end
%Create Yn
for k=1:K
    y(k)= x(k) + sqrt(sy)*randn();
end
%% RMSE
max_y = 0.1+max(y);
min_y = min(y)+0.2*(max_y-min(y));
th_s  = linspace(min_y,max_y,20);
save('compass_example_var','th_s','sx','sy','x','y')%save th_s, x, sx, sy

%% Set Model Parameters
In = ones(K,1);
Yn = y';
obs_valid=ones(K,1)*2;
save('model4_variables','In','Yn','obs_valid')

%% Create Original Model
compass_create_model

load('em_set100')
Param.Ak = 0.95;
save('em_set100')