%% INSTRUCTIONS:
%Running this script will call compass_create_model three times (once for 
%each model). Each time that this function is run you will be required to 
%answer prompted questions from the command window. In order to easily
%answer these questions, guided comments are provided which will take you
%through the entire compass_create_model function. These comments are
%located right below the line of code that compass_create_model is called,
%under each "Model #" heading.

%Once all three models are created, compass_run_models will be called and
%each created model will be run through compass_em for 250 iterations. It
%is important that you follow the guided comments exactly as they are
%written when answering questions prompted by compass_create_model, or else
%you will generate errors.






%% Model #1

%%Load behavioral data and prepare it for the toolbox
% Data: Yb , Yn
load('example_a.mat');
% Yn - logarithm of reaction time

ind = find(~isnan(Yn));               
ind2 = find(isnan(Yn));  %find nan values
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
obs_valid = ones(N,1);                            
obs_valid(ind2)=0; %set NaN values equal to zero 

%save variable set in current directory
save('model1_variables','Uk','In','Ib','Yn','Yb','obs_valid')


compass_create_model() %% Build behavioral model and learning procedure

%% The following comments are copied from compass_create_model input

% Please provide the exact path to your mat variable input file: model1_variables.mat
% 
%TIP (not from function): Merely entering the filename will only work if the file is in your
%%current directory.
%
% Do you have continuous observations? (Y/N)y
% 
% Do you have discrete observations? (Y/N)y
% 
% Which distribution for continuous observation?
% 
% Press 1 for normal distribution, or 2 for gamma distribution: 1
% 
% How many state variables do you have (nx = 1,2,3,...)1
% 
% *CONTINUOUS STATE VARIABLES*
% 
% Which state variable will be linked to column number 1 in your In matrix? (0=none,1,...,1) 1
% Is this a free or fixed parameter? (In(1)* b *x1) Enter 1 for free, 0 for fixed : 1
% What is the initial value of the parameter?1
% 
% Which state variable will be linked to column number 2 in your In matrix? (0=none,1,...,1) 0
% Is this a free or fixed parameter? (In(2)*b*) Enter 1 for free, 0 for fixed : 0
% What is the initial value of the parameter?0
% 
% *DISCRETE STATE VARIABLES*
% 
%  Which state variable will be linked to column number 1 of Ib ? (0=none,1,...,1) 1
% Is this a free or fixed parameter? (Ib(1)* b *x1) Enter 1 for free, 0 for fixed :0
% What is the initial value of the parameter? 1
% 
%  Which state variable will be linked to column number 2 of Ib ? (0=none,1,...,1) 0
% Is this a free or fixed parameter? (Ib(2)*b*) Enter 1 for free, 0 for fixed : 0
% What is the initial value of the parameter? 0
% 
% *The following 9 questions concern parameters for "compass_set_learning_param"* 
% 
% Enter "1" for "yes", or "0" for "no" *UNLESS OTHERWISE SPECIFIED*
% 
% Update state-transition model parameters?0
% Update state-transition model covariance matrix?1
% Update initial state variable parameters?0
% Update continuous parameters?1
% Update noise-term?1
% Update discrete parameters?1
% Make "A matrix" in state-transition process diagonal?1
% Update mean before covariance (ENTER 1), or covariance before mean (ENTER 2)? 1
% Update a positive shift estimation during the training process?0
% 
% Enter the number of compass_em training iterations.250
% Please enter the full path of where to save your new variable set, ending with a backslash (""): 
% Please enter the suffix you would like to add to the "em_set" filename: 1


%TIP (not from function): Save the mat file with a path and filename you can remember. If you do not
%specify a path (and simply push enter), the file will be saved in your current working directory

clear
%% Model #2

%This will be an identical model to #1, except it excludes discrete
%observation.

%%Load behavioral data and prepare it for the toolbox
% Data: Yb , Yn
load('example_a.mat');
% Yn - logarithm of reaction time

ind = find(~isnan(Yn)); %find nan values           
ind2 = find(isnan(Yn));  
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
obs_valid = ones(N,1);                            


obs_valid(ind2)=0; %set NaN values equal to zero 
%save variable set in current directory
save('model2_variables','Uk','In','Ib','Yn','Yb','obs_valid')
% %% Build behavioral model and learning procedure
compass_create_model()

%% The following comments are copied from compass_create_model input

% Please provide the exact path to your mat variable input file: model2_variables.mat
% 
%TIP (not from function): Merely entering the filename will only work if the file is in your
%%current directory.
%
% Do you have continuous observations? (Y/N)y
% 
% Do you have discrete observations? (Y/N)n
% 
% Which distribution for continuous observation?
% 
% Press 1 for normal distribution, or 2 for gamma distribution: 1
% 
% How many state variables do you have (nx = 1,2,3,...)1
% 
% *CONTINUOUS STATE VARIABLES*
% 
% Which state variable will be linked to column number 1 in your In matrix? (0=none,1,...,1) 1
% Is this a free or fixed parameter? (In(1)* b *x1) Enter 1 for free, 0 for fixed : 1
% What is the initial value of the parameter?1
% 
% Which state variable will be linked to column number 2 in your In matrix? (0=none,1,...,1) 0
% Is this a free or fixed parameter? (In(2)*b*) Enter 1 for free, 0 for fixed : 0
% What is the initial value of the parameter?0
% 
% *The following 9 questions concern parameters for "compass_set_learning_param"* 
% 
% Enter "1" for "yes", or "0" for "no" *UNLESS OTHERWISE SPECIFIED*
% 
% Update state-transition model parameters?0
% Update state-transition model covariance matrix?1
% Update initial state variable parameters?0
% Update continuous parameters?1
% Update noise-term?1
% Update discrete parameters?1
% Make "A matrix" in state-transition process diagonal?1
% Update mean before covariance (ENTER 1), or covariance before mean (ENTER 2)? 1
% Update a positive shift estimation during the training process?0
% 
% Enter the number of compass_em training iterations.250
% Please enter the full path of where to save your new variable set, ending with a backslash (""): 
% Please enter the suffix you would like to add to the "em_set" filename: 2

%TIP (not from function): Save the mat file with a path and filename you can remember. If you do not
%specify a path (and simply push enter), the file will be saved in your current working directory

clear
%% Model #3
%%Load behavioral data and prepare it for the toolbox
% Data: Yb , Yn
load('example_a.mat');
Yn = Yn(:,1);
Yb = Yb(:,1);
% Yn - logarithm of reaction time
ind = find(~isnan(Yn)); %find nan values               
ind2 = find(isnan(Yn)); 
ind3 = find(isnan(Yb));
Yn(ind2)=0;
Yo  = Yn/1000;
Yn = log(Yo);
% thresh=log(max(Yo)+.2);
Yn(ind2)=0;
Yb(ind3)=0;
N   = length(Yn);
% Input - 1 xi
In = zeros(N,2);In(:,1)= 1; In(:,2)= 1;
% Input, Ib is equal to In
Ib = In;
% Uk, which is all zero
Uk = zeros(N,1);
% valid, which is valid for the observed data point
obs_valid = ones(N,1);                            
obs_valid(ind2)=2;
In = In(:,1);
Ib = Ib(:,1);
%save variable set in current directory
save('model3_variables','Uk','In','Ib','Yn','Yb','obs_valid')

%%Build behavioral model and learning procedure
compass_create_model()

%% The following comments are copied from compass_create_model input

% Please provide the exact path to your mat variable input file: model3_variables.mat

%TIP (not from function): Merely entering the filename will only work if the file is in your
%%current directory.


% Do you have continuous observations? (Y/N)y
% 
% Do you have discrete observations? (Y/N)y
% 
% Which distribution for continuous observation?
% 
% Press 1 for normal distribution, or 2 for gamma distribution: 1
% 
% How many state variables do you have (nx = 1,2,3,...)1
% 
% *CONTINUOUS STATE VARIABLES*
% 
% Which state variable will be linked to column number 1 in your In matrix? (0=none,1,...,2) 1
% Is this a free or fixed parameter? (In(1)* b *x1) Enter 1 for free, 0 for fixed : 1
% What is the initial value of the parameter?1
% 
% CONTINUOUS MODEL: Yn~f(In(1)*(b1->1)*x1))
% Press enter to continue
% 
% *DISCRETE STATE VARIABLES*
% 
%  Which state variable will be linked to column number 1 of Ib ? (0=none,1,...,2) 1
% Is this a free or fixed parameter? (Ib(1)* b *x1) Enter 1 for free, 0 for fixed :0
% What is the initial value of the parameter? 1
%  
% DISCRETE MODEL: Yb~g(Ib(1)*1*x1))
% Press enter to continue
% 
% *The following 9 questions concern parameters for "compass_set_learning_param"* 
% 
% Enter "1" for "yes", or "0" for "no" *UNLESS OTHERWISE SPECIFIED*
% 
% Update state-transition model parameters?0
% Update state-transition model covariance matrix?1
% Update initial state variable parameters?0
% Update continuous parameters?1
% Update noise-term?1
% Update discrete parameters?1
% Make "A matrix" in state-transition process diagonal?1
% Update mean before covariance (ENTER 1), or covariance before mean (ENTER 2)? 1
% Update a positive shift estimation during the training process?0
% 
% Would you like to include a threshold in your model?(Y/N)y
% 
% Please enter a numerical value for your threshold. -.9
% 
% Please enter 1 to process censored data with imputation, 
% or enter 2 to process censored data with a Gaussian aproximation. 2
% 
% Enter 1 to update filter mean before covariance, or enter 2 to update covariance before filter mean.1
% 
% Enter the number of compass_em training iterations.250
% Please enter the full path of where to save your new variable set, ending with a backslash (""): 
% Please enter the suffix you would like to add to the "em_set" filename: 3

%TIP (not from function): Save the mat file with a path and filename you can remember. If you do not
%specify a path (and simply push enter), the file will be saved in your current working directory


%% compass_run_models


compass_run_models({'em_set1','em_set2','em_set3'})
%Entering the mat file names in this order will result in compass_em output
%saved (respectively) as model_result1, model_result2, and model_result3.








