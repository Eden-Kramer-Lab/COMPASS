function [Yn,Yb]=compass_sampling(DISTR,Cut_Time,Uk,In,Ib,Param,XPos0,SPos0)

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
% YP reaction time
% YB binary decision
if nargin  < 8    
    disp('insufficient input')
    return;
end

   
%% Build Mask Ck, Dk ,EK and Fk - note that Ck, Ek are time dependent and the Dk and Fk is linked to a subset of Input
Yn = [];
Yb = [];
if DISTR(1)>=1
    [MCk,MDk] = compass_Tk(In,Param);
    
    %% Normal Observation Model (  Y(k)=(Ck.*Tk)*X(k)+Dk*Ik+Vk*iid white noise )
    % ------------------
    % Y(k) is the observation, and Ik is the input either indicator or continuous
    % Ck, Dk, Vk are the model paramateres
    % Tk is model specific function - it is original set to but a one matrix
    % ------------------
    % Ck, NxM matrix - (Y is an observation of the length N, N can be 1, 2, ... - The Tk has the same size of input, 
    % and it is specfically designed for our task. It can be set to all 1 matrix)
    Ck = Param.Ck;           
    % Bk, NxS3 matrix - (We have an input of the length S3, and Dk will be size of NxS3) 
    Dk = Param.Dk.*MDk;           
    % Vk, is NxS4 matrix - (This is observation noise; for this, the S4 will be generally equal to N) 
    Vk = Param.Vk;
    % S
    if DISTR(1)==2
        S = Param.S;
    end
end
if DISTR(2)>=1
    [MEk,MFk] = compass_Qk(Ib,Param);

    %% Binary Observation Model (  P(k)=sigmoid((Ek.*Qk)*X(k)+Fk*Ik) )
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

%% State Space Model (  X(k+1)=Ak*X(k)+Bk*Uk+Wk*iid white noise )
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

% xMapping
xM = Param.xM;
    
%% Filtering Section
XPre = Ak * XPos0 + Bk * Uk';
SPre = Ak * SPos0* Ak' + Wk;

% Data observation: Normal
if DISTR(1) == 1   % main observation is Normal
     % Filtering
     CTk     = (Ck.*MCk{1})*xM;
     DTk     = Dk;
     % YP, Sk
     Sk      =  CTk * SPre * CTk' + Vk ;
     Yp      =  CTk * XPre + DTk * In';
     % Generate a sample - we assume it is scalar
     ys = Cut_Time:0.01:max(Cut_Time+10,Yp+10*sqrt(Sk)); 
     Pa  = pdf('normal',ys,Yp,sqrt(Sk));
     CPa = cumsum(Pa);
     CPa = CPa / sum(CPa);
     [~,ui] = min(abs(rand-CPa));
     Yn = ys(ui);
end
if DISTR(1) == 2   % main observation is Gamma
     % Filtering
     CTk     = (Ck.*MCk{1})*xM;
     DTk     = Dk;
     % YP, Sk
     EYn  =  exp(CTk * XPre + DTk * In');
     % Generate a sample - we assume it is scalar
     ys = Cut_Time:0.01:max(Cut_Time+10,EYn+5*EYn*EYn/Vk); 
     ys  = ys-S;
     Pa  = gampdf(ys,EYn*Vk,1/Vk);
     CPa = cumsum(Pa);
     CPa = CPa / sum(CPa);
     [~,ui] = min(abs(rand-CPa));
     Yn = S+ys(ui);
end
if DISTR(2)==1
    ETk  = (Ek.*MEk{1})*xM;
    FTk  = Fk;
    % calculate p
    temp = ETk * XPre + FTk * Ib';
    pk   = exp(temp)/(1+exp(temp));
    Yb   = 0;
    if rand() < pk 
        Yb=1;
    end
end
   
end
    
    