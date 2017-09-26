function [Rs,P]=compass_x_2_rp(DISTR,XP,Param,In,Ib)

%% Input Argument
    % MODE : leraning mode, defines how the learning value is updated
    %       MODE ==1, it means that the learning curve is determiend by X
    %       MODE ==2, it means that the learning curve is determiend by RT
    %       MODE ==3, it means that the learning curve is determiend by P - probability of true answer
    % Range: A vector definign upper and lower band values for the learning
    %       If it is based on X, it should be 2xlength(X), where the first row is the upper bound and the second is lowe bound (optional)
    %       If it is based on RT/P, it is a vector of 2x1, or scalar. Th first value defines the upper bound and the second is lowe band 
    % XP - Mean of the Estimate
    % SP - Covariance of the Estimate
    % Param, In, Ib - Check any other function
    
%% Output Argument
    % LRange: It is the learning rate either a vector of 2x1 or 1x1
    % The upper bound means what is the probability of being above upper bound
    % The lwer bound means what is the probability of being below lower bound


%% Reaction Time Observation - R
if DISTR(1)==1
    % Load Parameters
    Ck = Param.Ck;            
    Vk = Param.Vk;
    xM = Param.xM;
    [MCk,MDk] = compass_Tk(In,Param);
    DTk = Param.Dk.*MDk;           
    CTk = (Ck.*MCk{1})*xM;
    % Calculate Mean and Variance (we assume Log of RT)
    SSp  =  Vk ;
    MMp  =  CTk * XP + DTk * In';
    % Range 
    Rs(1)= icdf('Lognormal',0.975,MMp,sqrt(SSp));
    Rs(2)= icdf('Lognormal',0.025,MMp,sqrt(SSp));
    
end
%% Bonary Decision Observation
if DISTR(2)==1
    % Load Parameters
    Ek = Param.Ek;           
    [MEk,MFk] = compass_Qk(Ib,Param);
    FTk = Param.Fk.*MFk;           
    xM = Param.xM;
    ETk = (Ek.*MEk{1})*xM;
    % Calculate Mean and Variance 
    MMp=ETk * XP + FTk * Ib';
    P=exp(MMp)/(1+exp(MMp));
end    
    
    
end
    
    