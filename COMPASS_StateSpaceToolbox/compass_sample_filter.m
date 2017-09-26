function Xs =compass_sample_filter(Ns,DISTR,Uk,In,Ib,Yn,Yb,Param,obs_valid)
%% This draws samples from filter
K   = length(In);
d   = size(Param.Wk,1);
Xs  = zeros(Ns,K,d);
EPS = realmin('single');
MAX_EXP = 50;
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

Ak = Param.Ak;           
Bk = Param.Bk;           
Wk = Param.Wk;
xM = Param.xM;

%% Censored Reaction Time
if ~isempty(find(obs_valid==2))
    censor_time = Param.censor_time;
    Yn(find(obs_valid~=1))= censor_time;
end

%% Normal/Gamms Observation Model 
if DISTR(1)> 0
    Ck = Param.Ck;           
    Dk = Param.Dk.*MDk;           
    Vk = Param.Vk;
end

%% Binary Observation Model (P(k)=sigmoid((Ek.*Qk)*X(k)+Fk*Ik) )
if DISTR(2)==1
    Ek = Param.Ek;           
    Fk = Param.Fk.*MFk;           
end
%% Gamma, Extra Parameter - Time Shift
if DISTR(1)==2
    S = Param.S;
end
%% Check Uk
if isempty(Uk)
    Uk = zeros(K,size(Bk,2));
end

% Main Loop
for samp=1:Ns
    % display iter
    disp(['sample '  num2str(samp) ' out of ' num2str(Ns)])
    
    %% sample X0
    X0 = Param.X0 + chol(Param.W0)'*randn(d,1);  %mvnrnd(Param.X0,Param.W0)';
    W0 = 0*Param.W0;
    
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
                        st   = min(MAX_EXP,ETk * xpos + FTk * Ib(k,:)');
                        pk   = exp(st)./(1+exp(st));
                        xpos = XPre{k} +  SPre{k} * ETk' *(Yb(k)-pk);
                    end
                    XPos{k} = xpos;
                    % SPos
                    SPos{k} = (SPre{k}^-1 + ETk'*diag(pk.*(1-pk))*ETk)^-1;
                end
                % one-step mode
                if update_mode==2
                    st   =  min(MAX_EXP,ETk * XPre{k} + FTk * Ib(k,:)');
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
                            st   = min(MAX_EXP,ETk * xpos + FTk * Ib(k,:)');
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
                        st   =  min(MAX_EXP,ETk * XPre{k} + FTk * Ib(k,:)');
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
                            st   = min(MAX_EXP, ETk * xpos + FTk * Ib(k,:)');
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
                        st   = min(MAX_EXP, ETk * XPre{k} + FTk * Ib(k,:)');
                        pk   = exp(st)./(1+exp(st));
                        % SPos
                        SPos{k} = (SPre{k}^-1 + CTk' * CTk * Vk*(Yk/Yp) + ETk'*ETk*diag(pk.*(1-pk)) )^-1;
                        % XPos
                        XPos{k} = XPre{k} +  SPos{k} * ( ETk' *(Yb(k)-pk) - Vk * CTk'* (1-Yk/Yp) );
                    end
                end
            end
            %% Draw Samples
        else
            % randomly censored, the filter estimate will be equal to one step prediction
            XPos{k} = XPre{k};
            SPos{k} = SPre{k};
        end
        sXPos     = XPos{k} + chol(SPos{k})'*randn(d,1); %mvnrnd(XPos{k},SPos{k});
        XPos{k}   = sXPos;
        SPos{k}(:)= 0;
    end
    Xs(samp,:,:)=cell2mat(XPos')';
    
end
    
    