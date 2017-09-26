function [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb,MM]=compass_em_e(DISTR,Uk,In,Ib,Yn,Yb,Param,obs_valid)
%% Here, We assume each observation is scalar

%% Input Argument
    % DISTR, a vecotr of two variables. The [1 0] means there is only normal
    % observation/s, [0 1] means there is only binary observation/s, and [1 1]
    % will be both observations.
    % Uk: is a matrix of size KxS1 - K is the length of observation - input to
    % State model - X(k)=A*X(k-1)+B*Uk+Wk
    % In: is a matrix of size KxS3 - K is the length of observation - input to Normal observation model
    % Yn(k)=(C.*MCk)*X(k)+(D.*MDk)*In+Vk       - C and D are free parameters,
    % and MCk and MDk are input dependent components
    % Ib: is a matrix of size KxS5 - K is the length of observation - input to Binary observation model
    %  P(Yb(k)==1)=sigmoid((E.*MEk)*X(k)+(F.*MFk)*Ib       - E and F are free parameters,
    % and MEk and MFk are input dependent components
    % Yn: is a matrix of size KxN  - K is the length of observation, matrix of
    % normal observation
    % Yb: is a matrix of size KxN  - K is the length of observation, matrix of
    % binary observation
    % Param: it keeps the model information, and paramaters
%% Output Argument
    % XSmt is the smoothing result - mean
    % SSmt is the smoothing result - variance
    % Param is the updated model parameters
    % XPos is the filtering result - mean
    % SPos is the filtering result - variance
    % ML is the value of E-step maximization
    % EYn is the prdiction of the Yn
    % EYb is not added yet, but it can be the prediction of binary probability
%% obs_valid has three values (0= MAR, 1=observed, 2= censored)
EPS = realmin('single');

update_mode=Param.UpdateMode;
%% Observation Mode, from 1 to 5
if DISTR(1)==1
    observe_mode = DISTR(1) + 2*DISTR(2);
elseif DISTR(1)==2
    observe_mode = 2*DISTR(1) + DISTR(2);
else
    observe_mode = 2*DISTR(2);
end
%% Build Mask Ck, Dk ,EK and Fk - note that Ck, Ek are time dependent and the Dk and Fk is linked to a subset of Input
[MCk,MDk] = compass_Tk(In,Param);
if DISTR(2)==1
    [MEk,MFk] = compass_Qk(Ib,Param);
end

%% State Space Model (X(k+1)= Ak*X(k) + Bk*Uk + Wk*iid white noise )
% ------------------
% X(k) is the state, and Uk is the input
% Ak, Bk, Wk are model paramateres
% ------------------
% Ak, MxM matrix  (M is the length of the X)
Ak = Param.Ak;           
% Bk, MxS1 matrix (S1 is the length of Uk, Uk is a vector of size S1x1)
Bk = Param.Bk;           
% Wk, is MxS2 matrix (S2 is the length of Noise, we normally set the noise with the same dimension as the X - S2=M)
Wk = Param.Wk;
% X0, is a Mx1 matrix (initial value of X0)
X0 = Param.X0;
W0 = Param.W0;
% This is extending x
xM = Param.xM;

%% Censored Reaction Time
if ~isempty(find(obs_valid==2))
    censor_time = Param.censor_time;
    Yn(find(obs_valid~=1))= censor_time;
end

%% Normal/Gamms Observation Model 
if DISTR(1)> 0
    % For Normal,  Y(k)=(Ck.*Tk)*X(k)+Dk*Ik + Vk    Vk variance of iid white noise
    % For Gamma,   Y(k)=exp((Ck.*Tk)*X(k)+Dk*Ik)    Vk is dispersion term 
    % ------------------
    % Y(k) is the observation, and Ik is the input either indicator or continuous
    % Ck, Dk, Vk are the model paramateres
    % Tk is model specific function - it is original set to but a one matrix
    % ------------------
    % Ck, 1xM matrix - (Y is an scalar observation at each time point ... - The Tk has the same size of input, 
    % and it is specfically designed for our task. It can be set to all 1 matrix)
    Ck = Param.Ck;           
    % Bk, NxS3 matrix - (We have an input of the length S3, and Dk will be size of NxS3) 
    Dk = Param.Dk.*MDk;           
    % Vk, is scaler represnting noise in Normal or Dispresion Term in Gamma
    Vk = Param.Vk;
    
    % Length of data
    K    = length(Yn);
end

%% Binary Observation Model (P(k)=sigmoid((Ek.*Qk)*X(k)+Fk*Ik) )
if DISTR(2)==1
    % ------------------
    % P(k) is the observation probability at time k, and Ik is the input either indicator or continuous
    % Ek, and Fk are the model paramateres
    % Qk is model specific function - it is original set to but a one matrix
    % ------------------
    % Ck, NxM matrix - similar to Ck, Tk
    Ek = Param.Ek;           
    % Fk, NxS5 matrix - Similar to Dk
    Fk = Param.Fk.*MFk;           

    % Length of data
    K    = length(Yb);
end
%% Gamma, Extra Parameter - Time Shift
if DISTR(1)==2
    S = Param.S;
end

%% Check Uk
if isempty(Uk)
    Uk= zeros(K,size(Bk,2));
end

%% EM Update Loop
% Model Prediction
EYn = [];
EYb = [];
if DISTR(1)
   EYn  = zeros(K,1);
end
if DISTR(2)
   EYb  = zeros(K,1);
end
% ML array
ML = cell(Param.Iter,1);
% Main Loop
xx = -2.5:0.01:2.5;
aa = -2:0.05:0.5;
%aa = 0.01:0.01:1;
MM.CORRECT  = inf(Param.Iter,length(aa));
MM.CORMAX  = inf(Param.Iter,2);

MM.ESTIMATE = inf(Param.Iter,length(aa));
MM.ESTMAX = inf(Param.Iter,2);
MM.aa = aa;
%------
for iter=1:Param.Iter
    % display iter
    disp(['iteration '  num2str(iter) ' out of ' num2str(Param.Iter)])
    
    %% Run the Filtering Part
    % One step perdiction mean
    XPre = cell(K,1);
    % One step perdiction covariance
    SPre = cell(K,1);
    % Filter mean
    XPos = cell(K,1);
    % Filter covariance
    SPos = cell(K,1);
   
    % Filter 
    for k=1:K
        % One step prediction
        if k == 1
            XPre{k} = Ak * X0 + Bk * Uk(k,:)';
            SPre{k} = Ak * W0 * Ak'+ Wk;
        else
            XPre{k} = Ak * XPos{k-1} + Bk * Uk(k,:)';
            SPre{k} = Ak * SPos{k-1}* Ak' + Wk;
        end
        % Check if the data point is censored or not
        if obs_valid(k) 
            % Draw a sample if it is censored in censor_mode 1 (sampling)
            if obs_valid(k) == 2  && Param.censor_mode==1
                tIn  = [];
                tIb  = [];
                tUk  = [];
                if DISTR(1),        tIn = In(k,:); end
                if DISTR(2),        tIb = Ib(k,:); end
                if isempty(Uk)~=0,  tUk = Uk(k,:); end
                [tYP,tYB]=compass_sampling(DISTR,censor_time,tUk,tIn,tIb,Param,XPre{k},SPre{k});
                if DISTR(1),      Yn(k)=tYP;   end;
                if DISTR(2),      Yb(k)=tYB;   end;
            end
            % Observation: Normal
            if observe_mode == 1
                CTk     = (Ck.*MCk{k})*xM;
                DTk     =  Dk;
                if obs_valid(k) == 2  && Param.censor_mode==2
                    % censor time
                    T  = Param.censor_time;
                    if Param.censor_update_mode ==1
                        % SPos Update first
                        Mx = CTk * XPre{k} + DTk * In(k,:)';
                        Lx = max(EPS,normcdf(Yn(k),Mx,sqrt(Vk),'upper'));
                        Gx = normpdf(Yn(k),Mx,sqrt(Vk));
                        Tx = Gx/Lx;
                        % Calculate SPos
                        Hx = (Yn(k)-Mx)/Vk;
                        Sc = (CTk'*CTk)* Tx * (Tx-Hx);
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                        % XPos update next
                        Ac      =  CTk' * Tx;
                        XPos{k} =  XPre{k} + SPos{k} * Ac;
                    else
                        in_loop = 10;
                        % XPos update first
                        xpos = XPre{k};
                        for h= 1:in_loop
                            Mx = CTk * xpos + DTk * In(k,:)';
                            Lx = max(EPS,normcdf(Yn(k),Mx,sqrt(Vk),'upper'));
                            Gx = normpdf(Yn(k),Mx,sqrt(Vk));
                            Tx = Gx/Lx;
                            % update rule
                            Ac   = CTk' * Tx;
                            xpos = XPre{k} +  SPre{k} * Ac;
                        end
                        XPos{k} = xpos;
                        Mx = CTk * xpos + DTk * In(k,:)';
                        Lx = max(EPS,normcdf(Yn(k),Mx,sqrt(Vk),'upper'));
                        Gx = normpdf(Yn(k),Mx,sqrt(Vk));
                        Tx = Gx/Lx;
                        % SPos update next
                        Hx = (Yn(k)-Mx)/Vk;
                        Sc = (CTk'*CTk)*Tx * (Tx-Hx);
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                    end
                else
                    % XPos
                    Sk      =  CTk * SPre{k} * CTk' + Vk;
                    Yp      =  CTk * XPre{k} + DTk * In(k,:)';
                    XPos{k} =  XPre{k} + SPre{k} * CTk'* Sk^-1* (Yn(k)-Yp);
                    % SPos
                    SPos{k} = (SPre{k}^-1 + CTk' * Vk^-1 * CTk)^-1;
                end
            end
            % Observation: Bernoulli 
            % For Bernoulli mehtod, if there is any censored data it will be only based on resammpling technique
            if observe_mode == 2
                ETk = (Ek.*MEk{k})*xM;
                FTk = Fk;
                % XPos, SPos    
                % recursive mode
                if update_mode==1
                    in_loop = 10;
                    % XPos update
                    xpos = XPre{k};
                    for h= 1:in_loop
                        st   = ETk * xpos + FTk * Ib(k,:)';
                        pk   = exp(st)./(1+exp(st));
                        xpos = XPre{k} +  SPre{k} * ETk' *(Yb(k)-pk);
                    end
                    XPos{k} = xpos;
                    % SPos
                    SPos{k} = (SPre{k}^-1 + ETk'*diag(pk.*(1-pk))*ETk)^-1;
                end
                % one-step mode
                if update_mode==2
                    st   =  ETk * XPre{k} + FTk * Ib(k,:)';
                    pk   =  exp(st)./(1+exp(st));
                    SPos{k} = (SPre{k}^-1 + ETk'*diag(pk.*(1-pk))*ETk)^-1;
                    XPos{k} = XPre{k} +  SPos{k} * ETk' *(Yb(k)-pk);
                end
            end
            % Observation: Normal+Bernouli
            if observe_mode == 3
                CTk = (Ck.*MCk{k})*xM;
                DTk =  Dk;
                ETk = (Ek.*MEk{k})*xM;
                FTk =  Fk;
                if obs_valid(k) == 2  && Param.censor_mode==2
                    % This is exactly the same for Normal distribution
                    % censor time
                    T  = Param.censor_time;
                    % update mode 1
                    if Param.censor_update_mode==1
                        % SPos Update first
                        Mx = CTk * XPre{k} + DTk * In(k,:)';
                        Lx = max(EPS,normcdf(Yn(k),Mx,sqrt(Vk),'upper'));
                        Gx = normpdf(Yn(k),Mx,sqrt(Vk));
                        Tx = Gx/Lx;
                        % Update S
                        Hx = (Yn(k)-Mx)/Vk;
                        Sc = (CTk'*CTk)* Tx * (Tx-Hx);
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                        % XPos Update next
                        Ac      = CTk' * Tx;
                        XPos{k} =  XPre{k} + SPos{k} * Ac;
                    else
                        in_loop = 10;
                        % XPos update first
                        xpos = XPre{k};
                        for h= 1:in_loop
                            Mx = CTk * xpos + DTk * In(k,:)';
                            Lx = max(EPS,normcdf(Yn(k),Mx,sqrt(Vk),'upper'));
                            Gx = normpdf(Yn(k),Mx,sqrt(Vk));
                            Tx = Gx/Lx;
                            % S update
                            Ac      = CTk' * Tx;
                            xpos = XPre{k} +  SPre{k} * Ac;
                        end
                        XPos{k} = xpos;
                        Mx = CTk * xpos + DTk * In(k,:)';
                        Lx = max(EPS,normcdf(T,Mx,sqrt(Vk),'upper'));
                        Gx = normpdf(T,Mx,sqrt(Vk));
                        Tx = Gx/Lx;
                        % SPos update next
                        Hx = (Yn(k)-Mx)/Vk;
                        Sc = (CTk'*CTk)* Tx * (Tx-Hx);
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                    end
                else
                    % XPos, SPos
                    % recursive mode
                    if update_mode==1
                        % recursive mode
                        in_loop = 10;
                        xpos = XPre{k};
                        Yp   =  CTk * XPre{k} + DTk * In(k,:)';
                        Sk   =  (CTk' * Vk^-1 * CTk + SPre{k}^-1);
                        for z= 1:in_loop
                            st   = ETk * xpos + FTk * Ib(k,:)';
                            pk   = exp(st)./(1+exp(st));
                            xpos = XPre{k} +  Sk^-1 * ( ETk' *(Yb(k)-pk) + CTk'* Vk^-1 *(Yn(k)-Yp));
                        end
                        XPos{k} = xpos;
                        % SPos
                        SPos{k} = (SPre{k}^-1 + CTk' * Vk^-1 * CTk + ETk'*diag(pk.*(1-pk))*ETk )^-1;
                    end
                    % one-step mode
                    if update_mode==2
                        Yp   =  CTk * XPre{k} + DTk * In(k,:)';
                        st   =  ETk * XPre{k} + FTk * Ib(k,:)';
                        pk   =  exp(st)./(1+exp(st));
                        SPos{k} = (SPre{k}^-1 + ETk'*diag(pk.*(1-pk))*ETk + CTk' * Vk^-1 * CTk )^-1;
                        XPos{k} = XPre{k} +  SPos{k} * (ETk' *(Yb(k)-pk) + CTk'* (Yn(k)-Yp) * Vk^-1);
                    end
                end
            end
            % Observation: Gamma
            if observe_mode == 4
                CTk = (Ck.*MCk{k})*xM;
                DTk = Dk;
                % this is exactly equal to Normal case
                if obs_valid(k) == 2  && Param.censor_mode==2
                    % censor time 
                    if Param.censor_update_mode==1
                        % expected y
                        Mx = exp(CTk * XPre{k} + DTk * In(k,:)');
                        Hx = (Yn(k)-S)*Vk/Mx;
                        % components to estimate posterior
                        Lx = max(EPS,gammainc(Hx,Vk,'upper'));
                        Gx = gampdf(Hx,Vk,1);
                        % temporary
                        Ta = Gx/Lx;
                        % variace update
                        Sc = (CTk'*CTk)*((V-Hx)+Hx*Ta)*Hx*Ta;
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                        % XPos Update next
                        Ac      = CTk' *  Ta * Hx;
                        XPos{k} = XPre{k} + SPos{k} * Ac;
                    else
                        in_loop = 10;
                        % XPos update first
                        xpos = XPre{k};
                        for h= 1:in_loop
                            % expected y
                            Mx = exp(CTk * xpos + DTk * In(k,:)');
                            Hx = (Yn(k)-S)*Vk/Mx;
                            % components to estimate posterior
                            Lx = max(EPS,gammainc(Hx,Vk,'upper'));
                            Gx = gampdf(Hx,Vk,1);
                            % temporary
                            Ta   = Gx/Lx;
                            % XPos Update next
                            Ac   = CTk' *  Ta * Hx;
                            xpos = XPre{k} + SPre{k} * Ac;
                        end
                        XPos{k} = xpos;
                        Mx = exp(CTk * xpos + DTk * In(k,:)');
                        Hx = (Yn(k)-S)*Vk/Mx;
                        % components to estimate posterior
                        Lx = max(EPS,gammainc(Hx,Vk,'upper'));
                        Gx = gampdf(Hx,Vk,1);
                        % temporary
                        Ta = Gx/Lx;
                        % variace update
                        Sc = (CTk'*CTk)*((Vk-Hx)+Hx*Ta)*Hx*Ta;
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                    end
                else
                    % XPos, SPos
                    % recursive mode
                    if update_mode==1
                        % recursive mode
                        Yk      = Yn(k) - S;
                        in_loop = 10;
                        xpos    = XPre{k};
                        for h= 1:in_loop
                            Yp   = exp(CTk * xpos + DTk * In(k,:)');
                            xpos = XPre{k} - SPre{k} * Vk * CTk'* (1-Yk/Yp);    
                        end
                        XPos{k} = xpos;
                        SPos{k} = (SPre{k}^-1 + (Vk*(Yk/Yp))*CTk'*CTk)^-1;
                    end
                    if update_mode==2
                        Yk      = Yn(k) - S;
                        Yp      = exp(CTk * XPre{k} + DTk * In(k,:)');
                        SPos{k} = (SPre{k}^-1 + (Vk*(Yk/Yp))*CTk'*CTk)^-1;
                        XPos{k} =  XPre{k} - SPos{k} * Vk * CTk'* (1-Yk/Yp);
                    end
                end
            end
            % Observation: Gamma+Bernoulli
            if observe_mode == 5
                CTk = (Ck.*MCk{k})*xM;
                DTk = Dk;
                ETk = (Ek.*MEk{k})*xM;
                FTk = Fk;
                if obs_valid(k) == 2  && Param.censor_mode==2
                     % censor time 
                    if Param.censor_update_mode==1
                        % expected y
                        Mx = exp(CTk * XPre{k} + DTk * In(k,:)');
                        Hx = (Yn(k)-S)*Vk/Mx;
                        % components to estimate posterior
                        Lx = max(EPS,gammainc(Hx,Vk,'upper'));
                        Gx = gampdf(Hx,Vk,1);
                        % temporary
                        Ta = Gx/Lx;
                        % variace update
                        Sc = (CTk'*CTk)*((Vk-Hx)+Hx*Ta)*Hx*Ta;
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                        % XPos Update next
                        Ac      = CTk' *  Ta * Hx;
                        XPos{k} = XPre{k} + SPos{k} * Ac;
                    else
                        in_loop = 10;
                        % XPos update first
                        xpos = XPre{k};
                        for h= 1:in_loop
                            % expected y
                            Mx = exp(CTk * xpos + DTk * In(k,:)');
                            Hx = (Yn(k)-S)*Vk/Mx;
                            % components to estimate posterior
                            Lx = max(EPS,gammainc(Hx,Vk,'upper'));
                            Gx = gampdf(Hx,Vk,1);
                            % temporary
                            Ta = Gx/Lx;
                            % XPos Update next
                            Ac   = CTk' *  Ta * Hx;
                            xpos = XPre{k} + SPre{k} * Ac;
                        end
                        XPos{k} = xpos;
                        Mx = exp(CTk * xpos + DTk * In(k,:)');
                        Hx = (Yn(k)-S)*Vk/Mx;
                        % components to estimate posterior
                        Lx = max(EPS,gammainc(Hx,Vk,'upper'));
                        Gx = gampdf(Hx,Vk,1);
                        % temporary
                        Ta = Gx/Lx;
                        % variace update
                        Sc = (CTk'*CTk)*((Vk-Hx)+Hx*Ta)*Hx*Ta;
                        SPos{k} = ((SPre{k}^-1)+Sc)^-1;
                    end
                else
                    % recursive mode
                    if update_mode==1
                        % XPos, SPos
                        in_loop = 10;
                        Yk      = Yn(k) - S;
                        xpos    = XPre{k};
                        for h   = 1:in_loop
                            st   =  ETk * xpos + FTk * Ib(k,:)';
                            pk   = exp(st)./(1+exp(st));
                            Yp   = exp(CTk * xpos + DTk * In(k,:)');
                            xpos = XPre{k} +  SPre{k} * ( ETk' *(Yb(k)-pk) - Vk * CTk'* (1-Yk/Yp) );
                        end
                        XPos{k} = xpos;
                        % SPos
                        SPos{k} = (SPre{k}^-1 + CTk' * CTk * Vk *(Yk/Yp) + ETk'*ETk *diag(pk.*(1-pk)))^-1;
                    end
                    % one-step mode
                    if update_mode==2
                        % XPos, SPos
                        Yk   = Yn(k) - S;
                        Yp   = exp(CTk * XPre{k} + DTk * In(k,:)');
                        % Pk
                        st   =  ETk * XPre{k} + FTk * Ib(k,:)';
                        pk   = exp(st)./(1+exp(st));
                        % SPos
                        SPos{k} = (SPre{k}^-1 + CTk' * CTk * Vk*(Yk/Yp) + ETk'*ETk*diag(pk.*(1-pk)) )^-1;
                        % XPos
                        XPos{k} = XPre{k} +  SPos{k} * ( ETk' *(Yb(k)-pk) - Vk * CTk'* (1-Yk/Yp) );
                    end
                end
            end
        else
            % randomly censored, the filter estimate will be equal to one step prediction
            XPos{k} = XPre{k};
            SPos{k} = SPre{k};
        end
    end


    %% Smoother Part - it is based on classical Kalman Smoother
    % Kalman Smoothing
    As        = cell(K,1);
    % posterior mean
    XSmt      = cell(K,1);
    XSmt{end} = XPos{end};
    % posterior variance
    SSmt      = cell(K,1);
    SSmt{end} = SPos{end};
    for k=K-1:-1:1
        % Ak, equation (A.10)
        As{k} = SPos{k} * Ak' *  SPre{k+1}^-1 ;
        % Smting function, equation (A.9)
        XSmt{k} = XPos{k} + As{k} * (XSmt{k+1}- XPre{k+1});
        % Variance update, equation (A.11)
        SSmt{k} = SPos{k} + As{k} * (SSmt{k+1}- SPre{k+1}) * As{k}';
    end
    % Kalman smoother for time 0
    As0   = W0 * Ak' *  SPre{1}^-1;
    XSmt0 = X0 + As0 * (XSmt{1}- XPre{1}) ;
    SSmt0 = W0 + As0 * (SSmt{1}- SPre{1}) * As0';
    
    %% Extra Component of the State Prediction Ckk = E(Xk*Xk)
    % Ckk = E(Xk*Xk) prediction by smoothing 
    Ckk   = cell(K,1);
    for k = 1:K
        % Wk update - Smting Xk*Xk
        Ckk{k} = SSmt{k} + XSmt{k} * XSmt{k}';
    end
    Ckk0 = SSmt0 + XSmt0 * XSmt0';
    
    %% Extra Component of the State Prediction Ckk = E(Xk*Xk-1)
    % Ckk_1=E(Xk-1*Xk) prediction by smoothing - it is kept at index K
    % Wkk_1= Ckk_1 + Bias
    % Covariance for smoothed estimates in state space models - Biometrica
    % 1988- 601-602
    Ckk_1    = cell(K,1);
    Wkk_1    = cell(K,1);
    for k = 1:K
        % Wkk update - Smoothing Xk-1*Xk
        if k>1
            Wkk_1{k}   = As{k-1} * SSmt{k};
            Ckk_1{k}   = Wkk_1{k} + XSmt{k} * XSmt{k-1}'; 
        else
            Wkk_1{k}   = As0 * SSmt{k};
            Ckk_1{k}   = Wkk_1{k} + XSmt{k} * XSmt0'; 
        end
    end
    
    %% These are the function return value
    rXSmt = XSmt;
    rSSmt = SSmt;
    rXPos = XPos;
    rSPos = SPos;
    rYb   = Yb;
    rYn   = Yn;
    
    %% Here, We Generate EYn and Yb Prediction 
    %% Generate EYn
    if DISTR(1)>0
       for k=1:K
           % Filtering
           CTk     = (Ck.*MCk{k})*xM;
           DTk     = Dk;
           % EYn
           if DISTR(1)==1
                temp    = CTk * XSmt{k} + DTk * In(k,:)';
                EYn(k)  = temp';
           else
                temp    = CTk * XSmt{k} + DTk * In(k,:)';
                EYn(k)  = S + exp(temp);
           end
       end
    end
    if DISTR(2)==1
         for k=1:K
             % Filtering
             ETk     = (Ek.*MEk{k})*xM;
             FTk     = Fk;
             % YP
             temp    =  ETk * XSmt{k} + FTk * Ib(k,:)';
             EYb(k)  =  exp(temp')./(1+exp(temp'));
         end
    end
    
    
    %% Parameter Estimation Section
    %% Description of model parameters
    % Note that at index K of above variables, we keep:
    % 1. XPos(k) and SPos(k) are the optimum estimates of the Xk given 1:k
    % 2. XPre(k) and SPre(k) are the optimum estimates of the Xk give 1:k-1
    % 3. XSmt(k) and SSmt(k) are the optimum estimates of the Xk give 1:K
    % 4. Wk(k) is the E(Xk*Xk' given 1:K) observation
    % 5. Wkk(k)  is the E(Xk-1*Xk' given 1:K) observation
    % note that Wkk(1), E(X0*X1)=X0*E(X1) as X0 is a deterministic term

    % Calculate State Transition Model Parameters
    % Update Ak and Bk - we assume iid noise terms
     upAk = Ak;      
     upBk = Bk;
     upWk = Wk;
     % dimension of Xk 
     dx = size(Ak,2);
     % dimension of input
     di = size(Bk,2);
     if Param.UpdateStateParam == 1
         if Param.DiagonalA == 0
              % Update the Ak, Bk row by row
              % We define At once
              At = zeros(dx+di,dx+di);
              for k = 1:K
                 if k==1
                     % Calculate At
                     SecA = Ckk0;
                     SecB = XSmt0*Uk(k,:);
                     SecC = Uk(k,:)'*XSmt0';
                     SecD = Uk(k,:)'*Uk(k,:);
                     At   = At + [[SecA SecB]
                                  [SecC SecD]];
                 else
                     % Calculate At
                     SecA = Ckk{k-1};
                     SecB = XSmt{k-1}*Uk(k,:);
                     SecC = Uk(k,:)'*XSmt{k-1}';
                     SecD = Uk(k,:)'*Uk(k,:);
                     At   = At + [[SecA SecB]
                                  [SecC SecD]];
                 end
              end
              % We define Bt per row
              for d=1:dx
                  % Build At, Bt
                  Bt = zeros(dx+di,1);
                  for k = 1:K
                     % Calculate Bt
                     temp = Ckk_1{k};   SecA = temp(:,d);
                     temp = XSmt{k};    SecB = temp(d)* Uk(k,:)';
                     Bt   = Bt + [SecA;SecB];
                  end
                  % Calculate At and Bt for the d-th row
                  T  = pinv(At)*Bt;
                  upAk(d,:) = T(1:dx);
                  upBk(d,:) = T(dx+1:end);
              end
         end
         if Param.DiagonalA == 1
              % Update the Ak, Bk row by row
              for d=1:dx
                  % Build At, Bt
                  At = zeros(1+di,1+di);
                  Bt = zeros(1+di,1);
                  for k = 1:K
                      if k==1
                         % Calculate At
                         temp = Ckk0;
                         SecA = temp(d,d);
                         temp = XSmt0;
                         SecB = temp(d)*Uk(k,:);
                         SecC = Uk(k,:)'*temp(d);
                         SecD = Uk(k,:)'*Uk(k,:);
                         At   = At + [[SecA SecB]
                                      [SecC SecD]];
                      else
                         % Calculate At
                         temp = Ckk{k-1};
                         SecA = temp(d,d);
                         temp = XSmt{k-1};
                         SecB = temp(d)*Uk(k,:);
                         SecC = Uk(k,:)'*temp(d);
                         SecD = Uk(k,:)'*Uk(k,:);
                         At   = At + [[SecA SecB]
                                      [SecC SecD]];
                      end
                      % Calculate Bt
                      temp = Ckk_1{k};   SecA = temp(d,d)';
                      temp = XSmt{k};    SecB = temp(d)* Uk(k,:)';
                      Bt   = Bt + [SecA;SecB];
                  end
                  % Calculate At and Bt for the d-th row
                  T  = pinv(At)*Bt;
                  upAk(d,d) = T(1);
                  upBk(d,:) = T(2:end);
              end
         end
     end
     % Calculate Wk - we assume Wk is diagonal
     for d=1:dx
         upWk(d,d) = 0;
         for k = 1:K
              if k==1
                  % add E(Xk*Xk)
                  temp = Ckk{k};   upWk(d,d) = upWk(d,d) + temp(d,d);
                  % add E(Xk-1*Xk-1)*A^2
                  temp = Ckk0;     upWk(d,d) = upWk(d,d) + upAk(d,:) * temp * upAk(d,:)';
                  % add Bk*Uk^2
                                   upWk(d,d) = upWk(d,d) + (upBk(d,:) * Uk(k,:)')^2;
                  % add -2*A*E(Xk-1*Xk)                  
                  temp = Ckk_1{k}; upWk(d,d) = upWk(d,d) - upAk(d,:) * (temp(:,d)+ temp(d,:)');
                  % add -2*B*U*E(Xk)
                  temp = XSmt{k};  upWk(d,d) = upWk(d,d) - 2 * (upBk(d,:) *  Uk(k,:)') * temp(d);
                  % add 2 *B*U*A*E(Xk-1)
                  temp = XSmt0;    upWk(d,d) = upWk(d,d) + 2*(upBk(d,:) * Uk(k,:)')*(upAk(d,:)*temp);
              else
                  % add E(Xk*Xk)
                  temp = Ckk{k};   upWk(d,d) = upWk(d,d) + temp(d,d);
                  % add E(Xk-1*Xk-1)*A^2
                  temp = Ckk{k-1}; upWk(d,d) = upWk(d,d) + upAk(d,:) * temp * upAk(d,:)';
                  % add Bk*Uk^2
                                   upWk(d,d) = upWk(d,d) + (upBk(d,:) * Uk(k,:)')^2;
                  % add -2*A*E(Xk-1*Xk)                  
                  temp = Ckk_1{k}; upWk(d,d) = upWk(d,d) - upAk(d,:) * (temp(:,d)+ temp(d,:)');
                  % add -2*B*U*E(Xk)
                  temp = XSmt{k};  upWk(d,d) = upWk(d,d) - 2 * (upBk(d,:) *  Uk(k,:)') * temp(d);
                  % add 2 *B*U*A*E(Xk-1)
                  temp = XSmt{k-1};upWk(d,d) = upWk(d,d) + 2*(upBk(d,:) * Uk(k,:)')*(upAk(d,:)*temp);
              end
         end
         upWk(d,d) = upWk(d,d)/K;
         %---------------------------------
         % Update State parameters
         %----------------------------------
         Ak   = upAk;
         Bk   = upBk;
         if Param.UpdateStateNoise == 1
             Wk   = upWk;
         end
     end
     % Calculate the X0 parameters
     if Param.UpdateStateX0 == 1
        X0 = XSmt0;
        W0 = Wk;
     end
     % -----------------------------------------------
     % Calculate likelihood function (Hidden Variable)
     % Constant terms are excluded
     % -----------------------------------------------
     MaxH =  0;
     for d=1:size(Ak,2)
         %-- first variance 
         MaxH = MaxH -0.5*K*log(Wk(d,d));
         %-- other terms
         TempH = 0;
         for k = 1:K
              if k==1
                  % add E(Xk*Xk)
                  temp = Ckk{k};    TempH = TempH + temp(d,d);
                  % add E(Xk-1*Xk-1)*A^2
                  temp = Ckk0;      TempH = TempH + Ak(d,:) * temp * Ak(d,:)';
                  % add Bk*Uk^2
                                    TempH = TempH + (Bk(d,:) * Uk(k,:)')^2;
                  % add -2*A*E(Xk-1*Xk)                  
                  temp = Ckk_1{k};  TempH = TempH - Ak(d,:) * (temp(:,d)+ temp(d,:)');
                  % add -2*B*U*E(Xk)
                  temp = XSmt{k};   TempH = TempH - 2 * (Bk(d,:) *  Uk(k,:)') * temp(d);
                  % add 2 *B*U*A*E(Xk-1)
                  temp = XSmt0;     TempH = TempH + 2*(Bk(d,:) * Uk(k,:)')*(Ak(d,:)*temp);
              else
                  % add E(Xk*Xk)
                  temp = Ckk{k};    TempH = TempH + temp(d,d);
                  % add E(Xk-1*Xk-1)*A^2
                  temp = Ckk{k-1};	TempH = TempH + Ak(d,:) * temp * Ak(d,:)';
                  % add Bk*Uk^2
                                    TempH = TempH + (Bk(d,:) * Uk(k,:)')^2;
                  % add -2*A*E(Xk-1*Xk)                  
                  temp = Ckk_1{k};  TempH = TempH - Ak(d,:) * (temp(:,d)+ temp(d,:)');
                  % add -2*B*U*E(Xk)
                  temp = XSmt{k};   TempH = TempH - 2 * (Bk(d,:) *  Uk(k,:)') * temp(d);
                  % add 2 *B*U*A*E(Xk-1)
                  temp = XSmt{k-1}; TempH = TempH + 2*(Bk(d,:) * Uk(k,:)')*(Ak(d,:)*temp);
              end
         end
         MaxH = MaxH -0.5 * TempH / Wk(d,d);
     end
        
     %----------------------------------------------
     % Update the Observation Model Parameters
     %----------------------------------------------
     MaxO = 0;
     if observe_mode ==1 || observe_mode==3
        % replace unobserved points with censored threshold
        lb = zeros(K,1);
        for k = 1:K
            ctk = (Ck.*MCk{k})*xM;
            dtk =  Dk;
            dy  =  (Yn(k)-(ctk*XSmt{k}+dtk*In(k,:)'))/sqrt(Vk);
            sy  =  (ctk*SSmt{k}*ctk')/Vk;
           
            for h=1:3
                 temp = exp(dy)*exp(0.5*(1+lb(k))^2*sy);
                 lb(k)=- temp/(temp+exp(0.5*sy*lb(k)^2));
            end
        end
        % if we update all parameters of the model
        if Param.UpdateCModelParam == 1 && Param.UpdateCModelNoise == 1
            if ~isempty([find(Param.cLinkUpdate)  find(MDk)])
                % generate index and matrixes for optimization
                c_fill_ind = find(Param.cLinkUpdate);
                ck = Ck;
                d_fill_ind = find(MDk);
                dk = Dk;
                % initial parameters
                p0  = [Vk  ck(c_fill_ind)  dk(d_fill_ind)];
                % define bounds
                lower_bound=[EPS   -1e3*ones(1,length(p0)-1)];
                upper_bound=1e3*ones(1,length(p0));
                % call optimization function
                %options = optimoptions('fmincon','Display','off','DiffMaxChange',100,'MaxIter',1000,'Algorithm','interior-point');
                [p_opt,temp] = fminsearchbnd(@NormalParamCDV,p0,lower_bound,upper_bound);
                MaxO = -temp;
                disp('Normal');
                % put the estimates back to model
                Vk = p_opt(1);
                Ck(c_fill_ind)= p_opt(2:1+length(c_fill_ind));
                Dk(d_fill_ind) = p_opt(2+length(c_fill_ind):end);
            end
        end
        % if we only update Parameters of the model
        if Param.UpdateCModelParam == 1 && Param.UpdateCModelNoise == 0
            if ~isempty([find(Param.cLinkUpdate)  find(MDk)])
                % generate index and matrixes for optimization
                c_fill_ind = find(Param.cLinkUpdate);
                ck = Ck;
                d_fill_ind = find(MDk);
                dk = Dk;
                % initial parameters
                p0  = [ck(c_fill_ind)  dk(d_fill_ind)];
                % define bounds
                lower_bound=-1e3*ones(1,length(p0));
                upper_bound= 1e3*ones(1,length(p0));
                % call optimization function
                %options = optimoptions('fmincon','Display','off','DiffMaxChange',100,'MaxIter',1000,'Algorithm','interior-point');
                [p_opt,temp]   = fminsearchbnd(@NormalParamCD,p0,lower_bound,upper_bound);
                MaxO = -temp;
                disp('Normal');
                sv = Vk;
                
                [junk,ind]=min(abs(aa-p_opt(1)));
                MM.ESTMAX(iter,1) = -temp;
                MM.ESTMAX(iter,2) = aa(ind);
                
                for ll=1:length(aa)
                    MM.ESTIMATE(iter,ll)=-feval(@NormalParamCD,aa(ll));
                end
                
                for ind_aa=1:length(aa)
                    f=0;
                    for ind_kk=1:K
                        z = val_ind(ind_kk);
                        % param on time
                        ctk = (ck.*MCk{z})*xM;
                        dtk =  dk; dtk(2)=aa(ind_aa);
                        if obs_valid(z)==1
                            dy = Yn(z)-(ctk*xx+dtk*In(z,:)');
                            px = normpdf(xx,XSmt{z},sqrt(SSmt{z}));
                            px = px/sum(px);
                            f  = f - 0.5*log(sv)- (dy.^2*px')/(2*sv);
                        end
                        if obs_valid(z)==2
                            px = normpdf(xx,XSmt{z},sqrt(SSmt{z}));
                            px = px/sum(px);
                            f  = f +log(normcdf(Yn(z),ctk*xx+dtk*In(z,:)',sqrt(sv),'upper'))*px';
                        end
                    end
                    MM.CORRECT(iter,ind_aa)= f;
                end
                
                [max_val,t_ind]=max(MM.CORRECT(iter,:));
                MM.CORMAX(iter,1)= max_val;
                MM.CORMAX(iter,2)= aa(t_ind);
                % put the estimates back to model
                Ck(c_fill_ind)= p_opt(1:length(c_fill_ind));
                Dk(d_fill_ind) = p_opt(1+length(c_fill_ind):end);
            end
        end
        % if we only update Noise Term
        if Param.UpdateCModelParam == 0 && Param.UpdateCModelNoise == 1
                p0 = Vk;
                % define bounds
                lower_bound= EPS;
                upper_bound= 1e3;
                % call optimization function
                %options = optimoptions('fmincon','Display','off','DiffMaxChange',100,'MaxIter',1000,'Algorithm','interior-point');
                [p_opt,temp]   = fminsearchbnd(@NormalParamV,p0,lower_bound,upper_bound);
                MaxO = -temp;
                disp('Normal');
                
                ck = Ck;
                dk = Dk;
                
                
                [junk,ind]=min(abs(aa-p_opt(1)));
                MM.ESTMAX(iter,1) = -temp;
                MM.ESTMAX(iter,2) = aa(ind);
                
                for ll=1:length(aa)
                    MM.ESTIMATE(iter,ll)=-feval(@NormalParamV,aa(ll));
                end
                
                for ind_aa=1:length(aa)
                    f=0;
                    for ind_kk=1:K
                        z = val_ind(ind_kk);
                        % param on time
                        ctk = (ck.*MCk{z})*xM;
                        dtk =  dk;
                        sv  = aa(ind_aa);
                        if obs_valid(z)==1
                            dy = Yn(z)-(ctk*xx+dtk*In(z,:)');
                            px = normpdf(xx,XSmt{z},sqrt(SSmt{z}));
                            px = px/sum(px);
                            f  = f - 0.5*log(sv)- (dy.^2*px')/(2*sv);
                        end
                        if obs_valid(z)==2
                            px = normpdf(xx,XSmt{z},sqrt(SSmt{z}));
                            px = px/sum(px);
                            f  = f +log(normcdf(Yn(z),ctk*xx+dtk*In(z,:)',sqrt(sv),'upper'))*px';
                        end
                    end
                    MM.CORRECT(iter,ind_aa)= f;
                end
                
                [max_val,t_ind]=max(MM.CORRECT(iter,:));
                MM.CORMAX(iter,1)= max_val;
                MM.CORMAX(iter,2)= aa(t_ind);
                % put the estimates back to model
                Vk = p_opt;
        end
    end
    
    % Gamma observation
    if observe_mode ==4 || observe_mode==5
        % continuous parameters update
        if Param.UpdateCModelParam == 1
            if Param.UpdateCModelNoise && Param.UpdateCModelShift     % Full model update
                % keep learning param
                c_fill_ind = find(Param.cLinkUpdate);
                ck = Ck;  
                d_fill_ind = find(MDk==1);
                dk = Dk;
                % initiate parameters
                p0         = [Vk S ck(c_fill_ind) dk(d_fill_ind)];
                % lower and upper bound (disperssion Shift Ck Dk)
                lower_bound=[1      0               -1e3*ones(1,length(c_fill_ind))  -1e3*ones(1,length(d_fill_ind))];
                upper_bound=[50     0.99*min(Yn)          1e3*ones(1,length(c_fill_ind))   1e3*ones(1,length(d_fill_ind))];
                % call optimization function
                %options = optimoptions('lsqnonlin','Display','off','DiffMaxChange',100,'MaxIter',1000);
                [p_opt,temp]   = fminsearchbnd(@GammaParamFull,p0,lower_bound,upper_bound);
                MaxO = -temp;
                disp('Gamma');
                % update param
                Vk = p_opt(1);
                S  = p_opt(2); 
                Ck(c_fill_ind) = p_opt(3:2+length(c_fill_ind));
                Dk(d_fill_ind) = p_opt(3+length(c_fill_ind):end);
                
            elseif  Param.UpdateCModelNoise                           % Update Ck,Dk, plus V
                % keep learning param
                c_fill_ind = find(Param.cLinkUpdate);
                ck = Ck;  
                d_fill_ind = find(MDk==1);
                dk = Dk;
                % initiate p0
                p0         = [Vk ck(c_fill_ind) dk(d_fill_ind)];
                % lower and upper bound (disperssion Shift Ck Dk)
                lower_bound=[1      -1e3*ones(1,length(c_fill_ind))  -1e3*ones(1,length(d_fill_ind))];
                upper_bound=[50     1e3*ones(1,length(c_fill_ind))   1e3*ones(1,length(d_fill_ind))];
                % call optimization function
                %options = optimoptions('lsqnonlin','Display','off','DiffMaxChange',100,'MaxIter',1000);
                [p_opt,temp]   = fminsearchbnd(@GammaParamMinusS,p0,lower_bound,upper_bound);
                MaxO = -temp;
                disp('Gamma');
                
                % update param
                Vk = p_opt(1);
                Ck(c_fill_ind) = p_opt(2:1+length(c_fill_ind));
                Dk(d_fill_ind) = p_opt(2+length(c_fill_ind):end);
            elseif    Param.UpdateCModelShift                           % Update Ck,Dk, plus Shift
                 % keep learning param
                 c_fill_ind = find(Param.cLinkUpdate);
                 ck = Ck;  
                 d_fill_ind = find(MDk==1);
                 dk = Dk;
                 % initiate p0
                 p0         = [S ck(c_fill_ind) dk(d_fill_ind)];
                 % lower and upper bound (disperssion Shift Ck Dk)
                 lower_bound=[0               -1e3*ones(1,length(c_fill_ind))  -1e3*ones(1,length(d_fill_ind))];
                 upper_bound=[0.99*min(Yn)     1e3*ones(1,length(c_fill_ind))   1e3*ones(1,length(d_fill_ind))];
                 % call optimization function
                 %options = optimoptions('lsqnonlin','Display','off','DiffMaxChange',100,'MaxIter',1000);
                 [p_opt,temp]   = fminsearchbnd(@GammaParamMinusV,p0,lower_bound,upper_bound);
                 MaxO = -temp;
                 disp('Gamma');
                 % update param
                 S  = p_opt(1); 
                 Ck(c_fill_ind) = p_opt(2:1+length(c_fill_ind));
                 Dk(d_fill_ind) = p_opt(2+length(c_fill_ind):end);
            else
                % keep learning param
                c_fill_ind = find(Param.cLinkUpdate);
                ck = Ck;  
                d_fill_ind = find(MDk==1);
                dk = Dk;
                % initiate p0
                p0         = [ck(c_fill_ind) dk(d_fill_ind)];
                if ~isempty(p0)
                    % lower and upper bound (disperssion Shift Ck Dk)
                    lower_bound=[-1e3*ones(1,length(c_fill_ind))  -1e3*ones(1,length(d_fill_ind))];
                    upper_bound=[1e3*ones(1,length(c_fill_ind))   1e3*ones(1,length(d_fill_ind))];
                    % call optimization function
                    %options = optimoptions('lsqnonlin','Display','off','DiffMaxChange',100,'MaxIter',1000);
                    [p_opt,temp]   = fminsearchbnd(@GammaParamCD,p0,lower_bound,upper_bound);
                    MaxO = -temp;
                    disp('Gamma');
                    % update param
                    Ck(c_fill_ind) = p_opt(1:length(c_fill_ind));
                    Dk(d_fill_ind) = p_opt(1+length(c_fill_ind):end);
                end
            end
        end
    end
  
    %% Estimate Discrete Componentns
    MaxB = 0;
    if DISTR(2)==1
       if Param.UpdateDModelParam == 1
          %%-new update
          if ~isempty([find(Param.dLinkUpdate)  find(MFk)])
                % generate index and matrixes for optimization
                e_fill_ind = find(Param.dLinkUpdate);
                ek = Ek;
                f_fill_ind = find(MFk);
                fk = Fk;
                % initial parameters
                p0  = [ek(e_fill_ind)  fk(f_fill_ind)];
                % define bounds
                lower_bound=-1e3*ones(1,length(p0));
                upper_bound= 1e3*ones(1,length(p0));
                % calculate surrogate variables
                lb = zeros(K,1);
                for k = 1:K
                    etk = (ek.*MEk{k})*xM;
                    ftk =  fk;
                    ey  =  etk*XSmt{k}+ftk*Ib(k,:)';
                    sy  =  etk*SSmt{k}*etk';
                    for h=1:3
                        temp = exp(ey)*exp(0.5*(1+lb(k))^2*sy);
                        lb(k)=- temp/(temp+exp(0.5*sy*lb(k)^2));
                    end
                end
                % call optimization function
                %options = optimoptions('fmincon','Display','off','DiffMaxChange',100,'MaxIter',1000,'Algorithm','interior-point');
                [p_opt,temp] = fminsearchbnd(@BernoulliParam,p0,lower_bound,upper_bound);
                MaxB = -temp;
                disp('Bernoulli');
                % put the estimates back to model
                Ek(e_fill_ind) = p_opt(1:length(e_fill_ind));
                Fk(f_fill_ind) = p_opt(1+length(e_fill_ind):end);
          end
       end
    end
    ML{iter}.Total = MaxH + MaxO + MaxB;
    ML{iter}.State    = MaxH;
    ML{iter}.ObsrvNormGamma  = MaxO;
    ML{iter}.ObsrvBern  = MaxB;
end
%------------------------------
% Update Model Parameters
% State Transition Model Parameters
Param.Ak   = Ak;
Param.Bk   = Bk;
Param.Wk   = Wk;
Param.X0   = X0;
Param.W0   = W0;
% Binary section parameters
if DISTR(2) > 0
    Param.Ek   = Ek;
    Param.Fk   = Fk;
end
% Continiuous section parameters
if DISTR(1)> 0
    Param.Ck   = Ck;
    Param.Dk   = Dk;
    Param.Vk   = Vk;
    if DISTR(1)==2
        Param.S    = S;
    end
end
%-------------------------
    % Bernoulli Paramater Update Function
    function f= BernoulliParam(p)
        % replace param
        ek(e_fill_ind) = p(1:length(e_fill_ind));
        fk(f_fill_ind) = p(length(e_fill_ind)+1:end);
        % function calculation
        f=0;
        % now, calculate a part
        val_ind = find(obs_valid==1);
        for l = 1:length(val_ind)
           % valid indexes 
           z = val_ind(l);
           % param on time
           etk = (ek.*MEk{z})*xM;
           ftk =  fk;
           ey  =  etk*XSmt{z}+ftk*Ib(z,:)';
           sy  =  etk*SSmt{z}*etk';
           f   =  f - Yb(z)*ey + log(exp(0.5*lb(z)^2*sy)+exp(ey)*exp(0.5*(lb(z)+1)^2*sy));
        end
    end
    % Gamma Paramater Estimation
    function f = GammaParamFull(p)
        % first parameter is v
        % second parameter is alpha (alpha is positive and smaller than minimum Yn)
        % yk - paper EMBC
        yk = Yn - p(2);
        % define v
        v  = p(1);
        % ctk
        ck(c_fill_ind) = p(3:2+length(c_fill_ind));
        % dtk - parameters linked to Input
        dk(d_fill_ind) = p(3+length(c_fill_ind):end);
        % function value
        f = 0;
        % now, calculate a part
        val_ind = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           t = val_ind(l);
           % param on time
           ctk = (ck.*MCk{t})*xM;
           dtk = dk;
           % ey/sy
           ey  =  ctk*XSmt{t}+dtk*In(t,:)';
           sy  =  ctk*SSmt{t}*ctk';
           % commo term
           if obs_valid(t)==1
               f = f + log(gamma(v)) - v*log(yk(t)*v)+ v*ey + log(yk(t)) + v*yk(t)*exp(-ey+0.5*sy);
           end
           if obs_valid(t)== 2
               px = (yk(t)*exp(-ey)) / (1+yk(t)*exp(-ey));
               f  = f + v*yk(t)*exp(-ey+0.5*sy) - (v-1) * ( log(1+yk(t)*exp(-ey))+ 0.5*px*(1-px)*(ctk*SSmt{t}*ctk'));
           end
        end
    end
    % Gamma Paramater Estimation
    function f = GammaParamMinusV(p)
        % first parameter is v
        % second parameter is alpha (alpha is positive and smaller than minimum Yn)
        % yk - paper EMBC
        yk = Yn - p(1);
        % define v
        v  = Vk;
        % ctk
        ck(c_fill_ind) = p(2:1+length(c_fill_ind));
        % dtk - parameters linked to Input
        dk(d_fill_ind) = p(2+length(c_fill_ind):end);
        % function value
        f = 0;
        % now, calculate a part
        val_ind = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           t = val_ind(l);
           % param on time
           ctk = (ck.*MCk{t})*xM;
           dtk = dk;
           % ey/sy
           ey  =  ctk*XSmt{t}+dtk*In(t,:)';
           sy  =  ctk*SSmt{t}*ctk';
           % commo term
           if obs_valid(t)==1
               f = f + log(gamma(v)) - v*log(yk(t)*v)+ v*ey + log(yk(t)) + v*yk(t)*exp(-ey+0.5*sy);
           end
           if obs_valid(t)== 2
               px = (yk(t)*exp(-ey)) / (1+yk(t)*exp(-ey));
               f  = f + v*yk(t)*exp(-ey+0.5*sy) - (v-1) * ( log(1+yk(t)*exp(-ey))+ 0.5*px*(1-px)*(ctk*SSmt{t}*ctk'));
           end
        end
    end
    % Gamma Paramater Estimation
    function f = GammaParamMinusS(p)
        % first parameter is v
        % second parameter is alpha (alpha is positive and smaller than minimum Yn)
        % yk - paper EMBC
        yk = Yn - S;
        % define v
        v  = p(1);
        % ctk
        ck(c_fill_ind) = p(2:1+length(c_fill_ind));
        % dtk - parameters linked to Input
        dk(d_fill_ind) = p(2+length(c_fill_ind):end);
        % function value
        f = 0;
        % now, calculate a part
        val_ind = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           t = val_ind(l);
           % param on time
           ctk = (ck.*MCk{t})*xM;
           dtk = dk;
           % ey/sy
           ey  =  ctk*XSmt{t}+dtk*In(t,:)';
           sy  =  ctk*SSmt{t}*ctk';
           % commo term
           if obs_valid(t)==1
               f = f + log(gamma(v)) - v*log(yk(t)*v)+ v*ey + log(yk(t)) + v*yk(t)*exp(-ey+0.5*sy);
           end
           if obs_valid(t)== 2
               px = (yk(t)*exp(-ey)) / (1+yk(t)*exp(-ey));
               f  = f + v*yk(t)*exp(-ey+0.5*sy) - (v-1) * ( log(1+yk(t)*exp(-ey))+ 0.5*px*(1-px)*(ctk*SSmt{t}*ctk'));
           end
        end
    end
% Gamma Paramater Estimation
    function f = GammaParamCD(p)
        % first parameter is v
        % second parameter is alpha (alpha is positive and smaller than minimum Yn)
        % yk - paper EMBC
        yk = Yn - S;
        % define v
        v  = Vk;
        % ctk
        ck(c_fill_ind) = p(1:length(c_fill_ind));
        % dtk - parameters linked to Input
        dk(d_fill_ind) = p(1+length(c_fill_ind):end);
        % function value
        f = 0;
        % now, calculate a part
        val_ind = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           t = val_ind(l);
           % param on time
           ctk = (ck.*MCk{t})*xM;
           dtk = dk;
           % ey/sy
           ey  =  ctk*XSmt{t}+dtk*In(t,:)';
           sy  =  ctk*SSmt{t}*ctk';
           % commo term
           if obs_valid(t)==1
               f = f + log(gamma(v)) - v*log(yk(t)*v)+ v*ey + log(yk(t)) + v*yk(t)*exp(-ey+0.5*sy);
           end
           if obs_valid(t)== 2
               px = (yk(t)*exp(-ey)) / (1+yk(t)*exp(-ey));
               f  = f + v*yk(t)*exp(-ey+0.5*sy) - (v-1) * ( log(1+yk(t)*exp(-ey))+ 0.5*px*(1-px)*(ctk*SSmt{t}*ctk'));
           end
        end
    end

    % Normal parameters - C,D and Sv
    function f = NormalParamCDV(p)
        % logistic coef.
        la = -1.702;
        % replace params
        sv = p(1);
        ck(c_fill_ind)= p(2:1+length(c_fill_ind));
        dk(d_fill_ind)= p(2+length(c_fill_ind):end);
        % function value
        f  = 0;
        % now, calculate a part
        val_ind   = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           z = val_ind(l);
           % param on time
           ctk = (ck.*MCk{z})*xM;
           dtk =  dk;
           dy  =  Yn(z)-(ctk*XSmt{z}+dtk*In(z,:)');
           sy  =  ctk*SSmt{z}*ctk';
           if obs_valid(z)==1
                f = f + 0.5*log(sv)+ (dy^2+sy)/(2*sv);
           end
           if obs_valid(z)==2
                f = f - la * (dy/sqrt(sv)) + log( exp(0.5*lb(z)^2*sy/sv)+exp(la*dy/sqrt(sv))*exp(0.5*(lb(z)+la)^2*sy/sv));
           end
       end
    end
    % Normal parameters - C,D
    function f = NormalParamCD(p)
        % logistic coef.
        la = -1.702;
        % replace params
        sv = Vk;
        ck(c_fill_ind)= p(1:length(c_fill_ind));
        dk(d_fill_ind)= p(1+length(c_fill_ind):end);
        % function value
        f  = 0;
        % now, calculate a part
        val_ind   = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           z = val_ind(l);
           % param on time
           ctk = (ck.*MCk{z})*xM;
           dtk =  dk;
           dy  =  Yn(z)-(ctk*XSmt{z}+dtk*In(z,:)');
           sy  =  ctk*SSmt{z}*ctk';
           if obs_valid(z)==1
                f = f + 0.5*log(sv)+ (dy^2+sy)/(2*sv);
           end
           if obs_valid(z)==2
                f = f - la * (dy/sqrt(sv)) + log( exp(0.5*lb(z)^2*sy/sv)+exp(la*dy/sqrt(sv))*exp(0.5*(lb(z)+la)^2*sy/sv));
           end
       end
    end
    % Normal parameters - C,D
    function f = NormalParamV(p)
        % logistic coef.
        la = -1.702;
        % replace params
        ck = Ck;
        dk = Dk;
        sv = p(1);
        % function value
        f  = 0;
        % now, calculate a part
        % now, calculate a part
        val_ind   = find(obs_valid);
        for l =1:length(val_ind)
           % valid indexes 
           z = val_ind(l);
           % param on time
           ctk = (ck.*MCk{z})*xM;
           dtk =  dk;
           dy  =  Yn(z)-(ctk*XSmt{z}+dtk*In(z,:)');
           sy  =  ctk*SSmt{z}*ctk';
           if obs_valid(z)==1
                f = f + 0.5*log(sv)+ (dy^2+sy)/(2*sv);
           end
           if obs_valid(z)==2
                f = f - la * (dy/sqrt(sv)) + log( exp(0.5*lb(z)^2*sy/sv)+exp(la*dy/sqrt(sv))*exp(0.5*(lb(z)+la)^2*sy/sv));
           end
       end
    end
end
    
    