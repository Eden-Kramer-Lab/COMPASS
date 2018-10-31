function []= compass_create_model()
parameters = input('Please provide the exact path to your mat variable input file: ','s');
load(parameters)
%% Check that workspace variables exist
if (exist('In','var'))==0 || (exist('Yn','var'))==0
    fprintf('*WARNING: In and Yn do not concurrently exist* \nCheck uploaded data file to ensure necessary variables are included and spelled correctly.')
    continuous =0;
    In=[];
    Yn=[];
else
    continuous = 1;
end
if (exist('Ib','var'))==0||(exist('Yb','var'))==0
    fprintf('*WARNING: Ib and Yb do not concurrently exist* \nCheck uploaded data file to ensure necessary variables are included and spelled/capitalized correctly.')
    discrete=0;
    Ib=[];
    Yb=[];
else
    discrete = 1;
end
if (continuous == 0)&&(discrete==0)
    fprintf('*Error: Function cannot be executed due to lack of present variables* \n Check data file for necessary variables and proper spelling/capitalization.')
    return
end
if continuous
    K = size(In,1);
elseif discrete
    K = size(Ib,1);
end
if exist('obs_valid','var')==0
    obs_valid = ones(K,1);
end
if exist('Uk','var')==0
    Uk = [];
end
%% Initialize Variables
nUk=size(Uk,2);
nIn=size(In,2);
nIb=size(Ib,2);
thresh=9;
x=0;
censor_thr='g';
disc=[];
cont=[];
cont_mat=0;
disc_mat=0;
censor_mode       = 3;
update_mode=3;
UpdateStateParam  = 3;
UpdateStateNoise  = 3;
UpdateStateX0     = 3;
UpdateCModelParam = 3;
UpdateCModelNoise = 3;
UpdateDModelParam = 3;
DiagonalA         = 3;
UpdateMode        = 3;
UpdateCModelShift = 3;
%% Check for continuous/discrete observations
count=1; 
while ((strcmp(cont, 'n')==0)&&(strcmp(cont,'N')==0)&&(strcmp(cont,'y')==0)&&(strcmp(cont,'Y')==0))||((strcmp(disc, 'n')==0)&&(strcmp(disc,'N')==0)&&(strcmp(disc,'y')==0)&&(strcmp(disc,'Y')==0))
    if count>1
        disp('*Invalid answer to one or both questions*')
    end
    if continuous
    cont = input('\nDo you have continuous observations? (Y/N)','s');
    end
    if discrete
    disc = input('\nDo you have discrete observations? (Y/N)','s');
    end
    count=count+1;
end
if (strcmp(cont,'y')==1)||(strcmp(cont,'Y')==1)
    cont_mat=1;
    count=1;
    while ((x~=1)&&(x~=2))||(ischar(x)==1)
        if count>1
            fprintf('*Please enter either 1 or 2*')
        end
        x = input('\nWhich distribution for continuous observation?\n\nPress 1 for normal distribution, or 2 for gamma distribution: ');
        DISTR = [x 0];
        count=count+1;
    end
elseif (strcmp(cont, 'n')==1)||(strcmp(cont,'N')==1)
    DISTR = [0 0];
end
if ((disc == 'Y') || (disc =='y'))
    disc_mat=1;
    DISTR(2) = 1;
end
%% Determine state variables
nx    = input('\nHow many state variables do you have (nx = 1,2,3,...): ');
bc=0;
% link information
lIn  = zeros(nIn,1);
disc_str=[];
    % update information
    uIn  = zeros(nIn,1);
    % weight information
    wIn  = zeros(nIn,1);
if (cont_mat==1)

    % loop over n
    % Continuous state variable linking
    cont_str = [];
    for i=1:nIn %loop through In columns
        if i ==1
            fprintf('\n*CONTINUOUS STATE VARIABLES*\n')
        end
        count = 1;
        xv=99999999999999999;
        while xv > nx
            if count > 1
                fprintf('\n*Answer must be less than or equal to number of state variables present*')
            end
        str1 = ['\nWhich state variable will be linked to column number ' num2str(i) ' in your In matrix? (0=none,1,...,' num2str(nx) ') '];%link state variable to In column number
        xv = input(str1);
        count=count+1;
        end
        lIn(i)=xv;
        count=1;
        resp=9;
        while (resp~=1)&&(resp~=0)
            if count>1
                fprintf('*Please enter either 1 or 0*\n')
            end
        if xv > 0
            str1 = ['Is this a free or fixed parameter? (In(' num2str(i) ')*b*x'  num2str(xv) ') Enter 1 for free, 0 for fixed : '];
            
        else
            str1 = ['Is this a free or fixed parameter? (In(' num2str(i) ')*b*) Enter 1 for free, 0 for fixed : '];
            
        end
        resp = input(str1);
        count = count+1;
        end
        
        if resp ==1
            bc=bc+1;%beta count
        end
        uIn(i)=resp;
        str1 = 'What is the initial value of the parameter?';
        respo=input(str1);
        wIn(i)=respo;
        %Create Continuous Model String
        if (resp ==1) && (xv>0) %free parameter, with state variable
            cont_str=[cont_str '*(In(' num2str(i) ')*(b' num2str(bc)  '->' num2str(respo) ')*x'  num2str(xv) ')'];
        elseif (resp==0) && (xv>0) %fixed parameter, with state variable
            cont_str=[cont_str '*(In(' num2str(i) ')*' num2str(respo) '*x'  num2str(xv) ')'];
        elseif resp==1 %free parameter, without state variable
            cont_str=[cont_str '*(In(' num2str(i) ')*b' num2str(bc) '->' num2str(respo) ')'];
        else %fixed parameter, without state variable
            cont_str=[cont_str '*(In(' num2str(i) ')*' num2str(respo) ')'];
        end
    end
fprintf(['\nCONTINUOUS MODEL: Yn~f' cont_str(2:end) ')\n']) %display model
input('Press enter to continue')
end
wc=0; %initial w count
% link information
lIb  = zeros(nIb,1);
% update information
uIb  = zeros(nIb,1);
% weight informatio
wIb  = zeros(nIb,1);
if disc_mat==1
    % loop over n
    %Discrete observation variable linking
    for i=1:nIb%loop through Ib columns
        if i ==1
            fprintf('\n*DISCRETE STATE VARIABLES*\n')
        end
        count = 1;
        xv=9999999999999;
        while xv > nx
            if count > 1
                fprintf('\n*Answer must be less than or equal to number of state variables present*')
            end
        str1 = ['\n Which state variable will be linked to column number ' num2str(i) ' of Ib ? (0=none,1,...,' num2str(nx) ') '];%;link state variables to columns of In
        xv = input(str1);
        count=count+1;
        end
        lIb(i)=xv;
        count=1;
        resp=9;
        while (resp~=1)&&(resp~=0)
            if count>1
                fprintf('*Please enter either 1 or 0*\n')
            end
        if xv > 0
            str1 = ['Is this a free or fixed parameter? (Ib(' num2str(i) ')* w *x'  num2str(xv) ') Enter 1 for free, 0 for fixed : '];
        else
            str1 = ['Is this a free or fixed parameter? (Ib(' num2str(i) ')*w) Enter 1 for free, 0 for fixed : '];
        end
        resp = input(str1);
        count = count+1;
        end
        if resp ==1
            wc=wc+1;
        end
        uIb(i)=resp;
        str1 = 'What is the initial value of the parameter? ';
        respo=input(str1);
        wIb(i)=respo;
        if (resp ==1) && (xv>0) %free parameter, with state variable
            disc_str=[disc_str '*(Ib(' num2str(i) ')*(w' num2str(wc)  '->' num2str(respo) ')*x'  num2str(xv) ')'];
        elseif (resp==0) && (xv>0) %fixed parameter, with state variable
            disc_str=[disc_str '*(Ib(' num2str(i) ')*' num2str(respo) '*x'  num2str(xv) ')'];
        elseif resp==1 %free parameter, without state variable
            disc_str=[disc_str '*(Ib(' num2str(i) ')*w' num2str(wc) '->' num2str(respo) ')'];
        else %fixed parameter, without state variable
            disc_str=[disc_str '*(Ib(' num2str(i) ')*' num2str(respo) ')'];
        end
    end
    fprintf(['\nDISCRETE MODEL: Yb~g' disc_str(2:end) ')\n'])%display model
    input('Press enter to continue')
end








%Create xM matrix
xM = [];
for i=1:nx
    temp  = max(length(find(lIn==i)),length(find(lIb==i)));
    temp_a=  zeros(temp,nx); temp_a(:,i)=1;
    xM   = [xM;temp_a];
end
% cLink, cLinkUpdate, (Ck, DK) weight
cLink       = zeros(1,size(xM,1));
cLinkUpdate = zeros(1,size(xM,1));
Ck = ones(1,size(xM,1));
% SLink, dLinkUpdate, (Ek, Fk) weight
dLink       = zeros(1,size(xM,1));
dLinkUpdate = zeros(1,size(xM,1));
Ek = ones(1,size(xM,1));

txM = xM;
for i=1:length(lIn)
    if lIn(i)
        ind  = find(txM(:,lIn(i))==1);
        if isempty(ind)==0
            cLink(ind(1)) = i;
            txM(ind(1),:) = 0;
            cLinkUpdate(ind(1)) = uIn(i);
            Ck(ind(1)) = wIn(i);
        end
    end
end
ind = find(cLink==0);
cLink(ind)= nIn+(1:length(ind));
In  = [In zeros(size(In,1),length(ind))];
Dk  = zeros(1,size(In,2));
Dk(find(lIn==0)) = wIn(find(lIn==0))';
Dk(nIn+1:end)    = 1;
nIn = size(In,2);
txM = xM;
for  i=1:length(lIb)
    if lIb(i)
        ind  = find(txM(:,lIb(i))==1);
        if isempty(ind)==0
            dLink(ind(1)) = i;
            txM(ind(1),:) = 0;
            dLinkUpdate(ind(1)) = uIb(i);
            Ek(ind(1)) = wIb(i);
        end
    end
end
ind = find(dLink==0);
dLink(ind)= nIb+(1:length(ind));
Ib  = [Ib zeros(size(Ib,1),length(ind))];
Fk  = zeros(1,size(Ib,2));
Fk(find(lIb==0)) = wIb(find(lIb==0))';
Fk(nIb+1:end)    = 1;
nIb = size(Ib,2);
%% Run compass_create_state_space
Param    = compass_create_state_space(nx,nUk,nIn,nIb,xM,cLink,cLinkUpdate,dLink,dLinkUpdate);
Param.Ck = Ck;
Param.Dk = Dk;
Param.Ek = Ek;
Param.Fk = Fk;
%% Learning
learnstr = {'Update state-transition model parameters?' 'Update state-transition model covariance matrix?'...%Create cell of learning prompts
    'Update initial state variable parameters?' 'Update continuous parameters?' 'Update noise-term?' ...
    'Update discrete parameters?' 'Is matrix "A" in state-transition process diagonal or not?' ...
    'Update mean before covariance (ENTER 1), or covariance before mean (ENTER 2)? ' ...
    'Update a positive shift estimation during the training process?'};
fprintf('\n*Learning Parameter Update* \n\nEnter "1" for "yes", or "0" for "no" *UNLESS OTHERWISE SPECIFIED*\n\n')
a=3;
for i = 1:9 %loop through learning prompts
    count=1;
    if i ~= 8
        while (a~=1)&&(a~=0)
            if count>1
                disp('Only binary input is acceptable for this question')
            end
            a = input(char(learnstr(i)));
            if i ==1 %Each if and elseif statement corresponds to a prompt in the above cell
                UpdateStateParam = a;
            elseif i ==2
                UpdateStateNoise = a;
            elseif i ==3
                UpdateStateX0 = a;
            elseif (i ==4) && (cont_mat==1)
                UpdateCModelParam = a;
            elseif (i ==5) && (cont_mat==1)
                UpdateCModelNoise = a;
            elseif (i == 6) && (disc_mat==1)
                UpdateDModelParam = a;
            elseif i == 7
                DiagonalA=a;
            else
                UpdateCModelShift=a;
            end
            count=count+1;
        end
    else
        count = 1;
        while (a~=1)&&(a~=2)
            if count>1
                disp('Only 1 or 2 is a valid answer to this question')
            end
            a = input(char(learnstr(i)));
            count=count+1;
        end
        UpdateMode=a;
    end
    a=3;
end
fprintf('\nLearning parameters are set.\n')
input('Press "enter" to continue.')
%% Censoring Data Points
count=1;
if ismember(2,obs_valid)
    fprintf('\n*Threshold Parameters Update*\n')
while (thresh ~='y')&&(thresh~='Y')&&(thresh~='n')&&(thresh~='N')
    thresh=input('\nWould you like to include a threshold in your model?(Y/N)','s'); %include/exclude threshold
    if count>1
        fprintf('*Please enter a valid answer*')
    end
    count=count+1;
    if (thresh =='y')||(thresh=='Y')
        inc_thresh=1;
        while (ischar(censor_thr)==1)
            censor_thr = input('\nPlease enter a numerical value for your threshold: ');%enter threshold level
        end
        count=1;
        while (censor_mode~=1)&&(censor_mode~=2)
            if count>1
                disp('*Please enter a valid value: 1 or 2*')%choose censored data processing method
            end
            censor_mode = input('\nPlease enter 1 to process censored data with imputation, \nor enter 2 to process censored data with a Gaussian aproximation. ');
            count=count+1;
        end
    elseif (thresh=='n')||(thresh=='N')
        inc_thresh=0;
    else
    end
end
else
    inc_thresh=0;
end
if (thresh =='y')||(thresh=='Y')
if censor_mode==2
    count=1;
    while (update_mode~=1)&&(update_mode~=2)
        if count>1
            disp('*Please enter a valid value: 1 or 2*')
        end
        update_mode = input('\nEnter 1 to update filter mean before covariance, or enter 2 to update covariance before filter mean.'); %choose update_mode speed
        count=count+1;
    end
end
end
if ismember(2,obs_valid)
fprintf('\nThreshold parameters are set.\n')
input('Press enter to continue.')
end
%Run threshold
if inc_thresh
    Param = compass_set_censor_threshold_proc_mode(Param, censor_thr,censor_mode,update_mode);
end
Iter= input('\nEnter the number of compass_em training iterations: '); %Choose Iter value for compass_em
Param.Iter = Iter;
Param    = compass_set_learning_param(Param,Iter,UpdateStateParam,...
    UpdateStateNoise,UpdateStateX0,UpdateCModelParam,UpdateCModelNoise,...
    UpdateDModelParam,DiagonalA,UpdateMode,UpdateCModelShift);
path = input('Please enter the full path of where to save your new variable set, ending with a backslash ("\") \n(Press enter to save in current directory)','s');
suffix = input('Please enter the suffix you would like to add to the "em_set" filename: ','s');
name = ['em_set' suffix];
%% Save new variables
save([path name],'DISTR','Uk','In','Ib','Yn','Yb','Param','obs_valid');
end






