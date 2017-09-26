function [DEV_C,DEV_D]=compass_deviance(DISTR,In,Ib,Yn,Yb,Param,obs_valid,XSmt,SSmt)
% This function calculates the deviance for both continuous and discrete observations
% number of samples over censored points
Ns  = 10;
K   = max(length(Yn),length(Yb));
%% Here, you need the function to replace the censored/missing data with a sampleed one
if ~isempty(Yn)
    tYn = repmat(Yn,1,Ns); 
end
if ~isempty(Yb)
    tYb = repmat(Yb,1,Ns);
end

for k=1:K
    if obs_valid(k)==0 % missing 
        for s=1:Ns
            [samp_yn,samp_yb] = compass_post_sampling(DISTR,[],In(k,:),Ib(k,:),Param,XSmt{k},SSmt{k});
            if DISTR(1),   tYn(k,s)= samp_yn;   end
            if DISTR(2),   tYb(k,s)= samp_yb;   end
        end
    end
    if obs_valid(k)==2 % censored
        for s=1:Ns
            [samp_yn,samp_yb]   = compass_post_sampling(DISTR,Param.censor_time,In(k,:),Ib(k,:),Param,XSmt{k},SSmt{k});
            if DISTR(1),   tYn(k,s)= samp_yn;   end
            if DISTR(2),   tYb(k,s)= samp_yb;   end
        end
    end
end


%% CONTINUOUS PART DEVIANCE
% Normal
DEV_C = [];
if DISTR(1) == 1
    xM        = Param.xM;
    [MCk,MDk] = compass_Tk(In,Param);
    Vk   = Param.Vk;
    DTk  = Param.Dk.*MDk;           
    Ck   = Param.Ck; 
    % deviance calculation
    N     = 1000;
    DEV_C = 0;
    for  k=1:K
       % draw samples
       Xs  = mvnrnd(XSmt{k},SSmt{k},N)' ;  
       CTk = (Ck.*MCk{k})*xM;
       Mx = CTk * Xs + DTk * In(k,:)';
       Sx = Vk; 
       if obs_valid(k)==1 
            DEV_C = DEV_C -2*sum(log(pdf('normal',tYn(k,1),Mx,sqrt(Sx))))/N;
       end
       if obs_valid(k)==2 
            avg_log_ll = 0;
            for s=1:Ns
                avg_log_ll = avg_log_ll -2*sum(log(pdf('normal',tYn(k,s),Mx,sqrt(Sx))))/N;
            end
            avg_log_ll = avg_log_ll/Ns;
            DEV_C = DEV_C + avg_log_ll;
       end
    end
end
% Gamma
if DISTR(1)==2
    xM        = Param.xM;
    [MCk,MDk] = compass_Tk(In,Param);
    DTk = Param.Dk.*MDk;           
    % model parameters
    Ck = Param.Ck;
    Vk = Param.Vk;
    S  = Param.S;
    % draw samples per trial
    N    = 1000;
    DEV_C = 0;
    for  k=1:K
       Xs  = mvnrnd(XSmt{k},SSmt{k},N)' ; 
       CTk = (Ck.*MCk{k})*xM;
       Mx  = exp(CTk * Xs + DTk * In(k,:)');
       if obs_valid(k)==1
            DEV_C = DEV_C -2*sum(log(gampdf(tYn(k,1)-S,Vk,Mx/Vk)))/N;
       end
       if obs_valid(k)==2    
            avg_log_ll = 0;
            for s=1:Ns
                avg_log_ll = avg_log_ll -2*sum(log(gampdf(tYn(k,s)-S,Vk,Mx/Vk)))/N;
            end
            avg_log_ll = avg_log_ll/Ns;
            DEV_C = DEV_C + avg_log_ll;
       end
    end
end
%% Discrete Part DEVIANCE
% we assume fully observed data
DEV_D = [];
if DISTR(2)==1
    xM        = Param.xM;
    [MEk,MFk] = compass_Qk(Ib,Param);
    % model parameters
    Ek = Param.Ek;
    FTk = Param.Fk.*MFk;
    % draw samples per trial
    N    = 1000;
    % map to larger space
    DEV_D = 0;
    for k=1:K
           Xs  = mvnrnd(XSmt{k},SSmt{k},N)' ; 
           ETk = (Ek.*MEk{k})*xM;
           st = ETk * Xs + FTk * Ib(k,:)' ;
           pk = exp(st)./(1+exp(st));
           if obs_valid(k)== 1
                if tYb(k,1)
                     DEV_D = DEV_D -2*sum(log(pk))/N;
                else
                     DEV_D = DEV_D -2*sum(log(1-pk))/N;
                end
           end
           if obs_valid(k)== 2
               avg_log_ll = 0;
               for s=1:Ns
                 if tYb(k,s)
                     avg_log_ll = avg_log_ll -2*sum(log(pk))/N;
                 else
                     avg_log_ll = avg_log_ll -2*sum(log(1-pk))/N;
                 end
               end
               avg_log_ll = avg_log_ll/Ns;
               DEV_D = DEV_D + avg_log_ll;
          end
    end
end
    

