function [XPos,SPos,YP,YB]=compass_filtering(DISTR,Uk,In,Ib,Yn,Yb,Param,obs_valid,XPos0,SPos0)

%% Input Argument
    % DISTR, a vecotr of two variables. The [1 0] means there is only normal
    % observation/s, [0 1] means there is only binary observation/s, and [1 1]
    % will be both observations.
    % Uk: is a matrix of size KxS1 - K is the length of observation - input to
    % State model - X=A*X(k-1)+B*Uk+Wk
    % In: is a matrix of size KxS3 - K is the length of observation - input to Normal observation model
    % Yn=(C.*MCk)*X(k)+(D.*MDk)*In+Vk       - C and D are free parameters,
    % and MCk and MDk are input dependent components
    % Ib: is a matrix of size KxS5 - K is the length of observation - input to Binary observation model
    %  P(Yb==1)=sigmoid((E.*MEk)*X(k)+(F.*MFk)*Ib       - E and F are free parameters,
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
    % YP is the prdiction of the Yn
    % YB is not added yet, but it can be the prediction of binary probability
    
    
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
% This is extending x
xM = Param.xM;

%% Censored Reaction Time
if ~isempty(find(obs_valid==2))
    censor_time = Param.censor_time;
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

end

%% Check Uk
if isempty(Uk)
    Uk= zeros(1,size(Bk,2));
end

    
% Filter 

% One step prediction

XPre = Ak * XPos0 + Bk * Uk';
SPre = Ak * SPos0* Ak' + Wk;

% Check if the data point is censored or not
if ( obs_valid>=1 )
    % Draw a sample if it is censored
    if obs_valid==2    
        [tYP,tYB]=compass_sampling(DISTR,censor_time,Uk,In,Ib,Param,XPre,SPre);
        if DISTR(1)>0       Yn=tYP;   end;
        if DISTR(2)==1      Yb=tYB;   end;
    end
    % Observation: Normal
    if observe_mode == 1
        CTk     = (Ck.*MCk{1})*xM;
        DTk     =  Dk;
        % XPos
        Sk      =  CTk * SPre * CTk' + Vk;
        Yp      =  CTk * XPre + DTk * In';
        XPos =  XPre + SPre * CTk'* Sk^-1* (Yn-Yp);
        % SPos
        SPos = (SPre^-1 + CTk' * Vk^-1 * CTk)^-1;
    end
    % Observation: Bernoulli
    if observe_mode == 2
        ETk = (Ek.*MEk)*xM;
        FTk = Fk;
        % XPos, SPos    
        % recursive mode
        if update_mode==1
            in_loop = 10;
            % XPos update
            xpos = XPre;
            for t= 1:in_loop
                st   = ETk * xpos + FTk * Ib';
                pk   = exp(st)./(1+exp(st));
                xpos = XPre +  SPre * ETk' *(Yb-pk);
            end
            XPos = xpos;
            % SPos
            SPos = (SPre^-1 + ETk'*diag(pk.*(1-pk))*ETk)^-1;
        end
        % one-step mode
        if update_mode==2
            st   =  ETk * XPre + FTk * Ib';
            pk   =  exp(st)./(1+exp(st));
            SPos = (SPre^-1 + ETk'*diag(pk.*(1-pk))*ETk)^-1;
            XPos = XPre +  SPos * ETk' *(Yb-pk);
        end
    end
    % Observation: Normal+Bernouli
    if observe_mode == 3
        CTk = (Ck.*MCk{1})*xM;
        DTk =  Dk;
        ETk = (Ek.*MEk{1})*xM;
        FTk =  Fk;
        % XPos, SPos
        % recursive mode
        if update_mode==1
            % recursive mode
            in_loop = 10;
            xpos = XPre;
            Yp   =  CTk * XPre + DTk * In';
            Sk   =  (CTk' * Vk^-1 * CTk + SPre^-1);
            for t= 1:in_loop
                st   = ETk * xpos + FTk * Ib';
                pk   = exp(st)./(1+exp(st));
                xpos = XPre +  Sk^-1 * ( ETk' *(Yb-pk) + CTk'* Vk^-1 *(Yn-Yp));
            end
            XPos = xpos;
            % SPos
            SPos = (SPre^-1 + CTk' * Vk^-1 * CTk + ETk'*diag(pk.*(1-pk))*ETk )^-1;
        end
        % one-step mode
        if update_mode==2
            Yp   =  CTk * XPre + DTk * In';
            st   =  ETk * XPre + FTk * Ib';
            pk   =  exp(st)./(1+exp(st));
            SPos = (SPre^-1 + ETk'*diag(pk.*(1-pk))*ETk + CTk' * Vk^-1 * CTk )^-1;
            XPos = XPre +  SPos * (ETk' *(Yb-pk) + CTk'* (Yn-Yp) * Vk^-1);
        end
    end
    % Observation: Gamma
    if observe_mode == 4
        CTk = (Ck.*MCk{1})*xM;
        DTk = Dk;
        % XPos, SPos
        % recursive mode
        if update_mode==1
            % recursive mode
            Yk      = Yn - Param.S;
            in_loop = 10;
            xpos    = XPre;
            for t= 1:in_loop
                Yp   = exp(CTk * xpos + DTk * In');
                xpos = XPre - SPre * Vk * CTk'* (1-Yk/Yp);    
            end
            XPos = xpos;
            SPos = (SPre^-1 + (Vk*(Yk/Yp))*CTk'*CTk)^-1;
        end
        if update_mode==2
            Yk      = Yn - Param.S;
            Yp      = exp(CTk * XPre + DTk * In');
            SPos = (SPre^-1 + (Vk*(Yk/Yp))*CTk'*CTk)^-1;
            XPos =  XPre - SPos * Vk * CTk'* (1-Yk/Yp);
        end
    end
    % Observation: Gamma+Bernoulli
    if observe_mode == 5
        CTk = (Ck.*MCk{1})*xM;
        DTk = Dk;
        ETk = (Ek.*MEk{1})*xM;
        FTk = Fk;
        % recursive mode
        if update_mode==1
            % XPos, SPos
            in_loop = 10;
            Yk      = Yn - Param.S;
            xpos    = XPre;
            for t   = 1:in_loop
                st   =  ETk * xpos + FTk * Ib';
                pk   = exp(st)./(1+exp(st));
                Yp   = exp(CTk * xpos + DTk * In');
                xpos = XPre +  SPre * ( ETk' *(Yb-pk) - Vk * CTk'* (1-Yk/Yp) );
            end
            XPos = xpos;
            % SPos
            SPos = (SPre^-1 + CTk' * CTk * Vk *(Yk/Yp) + ETk'*ETk *diag(pk.*(1-pk)))^-1;
        end
        % one-step mode
        if update_mode==2
            % XPos, SPos
            Yk   = Yn - Param.S;
            Yp   = exp(CTk * XPre + DTk * In');
            % Pk
            st   =  ETk * XPre + FTk * Ib';
            pk   = exp(st)./(1+exp(st));
            % SPos
            SPos = (SPre^-1 + CTk' * CTk * Vk*(Yk/Yp) + ETk'*ETk*diag(pk.*(1-pk)) )^-1;
            % XPos
            XPos = XPre +  SPos * ( ETk' *(Yb-pk) - Vk * CTk'* (1-Yk/Yp) );
        end
    end
else
    % randomly censored, the filter estimate will be equal to one step prediction
    XPos = XPre;
    SPos = SPre;
end



%% Update prediction
YP = [];
YB = [];
if DISTR(1)>0
   
      % Filtering
      CTk     = (Ck.*MCk{1})*xM;
      DTk     = Dk;
      % EYn
      if DISTR(1)==1
          temp    = CTk * XPos + DTk * In';
          YP      = temp';
      else
          temp    = CTk * XPos + DTk * In';
          YP      = exp(0.5* CTk * temp * CTk');
      end

end
if DISTR(2)==1
   % Filtering
   ETk     = (Ek.*MEk{1})*xM;
   FTk     = Fk;
   % YP
   temp    =  ETk * XPos + FTk * Ib';
   YB  =  exp(temp')./(1+exp(temp'));
    
end


