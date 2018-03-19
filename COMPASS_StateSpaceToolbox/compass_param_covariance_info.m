function [COV_X,COV_C,COV_D]=compass_param_covariance_info(DISTR,Uk,In,Ib,Yn,Yb,Param,obs_valid,XSmt,SSmt)
% This function calculates covariance estimates of the model parameters

K   = max(length(Yn),length(Yb));
%% Here, you need the function to replace the censored/missing data with a sampleed one
for k=1:K
    if obs_valid(k)==0 % missing 
        [samp_yn,samp_yb] = compass_post_sampling(DISTR,[],In(k,:),Ib(k,:),Param,XSmt{k},SSmt{k});
        if DISTR(1),   Yn(k)=samp_yn;   end
        if DISTR(2),   Yb(k)=samp_yb;   end
    end
    if obs_valid(k)==2 % censored
        [samp_yn,samp_yb] = compass_post_sampling(DISTR,Param.censor_time,In(k,:),Ib(k,:),Param,XSmt{k},SSmt{k});
        if DISTR(1),   Yn(k)=samp_yn;   end
        if DISTR(2),   Yb(k)=samp_yb;   end
    end
end

%% CA, CB, CAB, CW - note we assume W is diagonal
for r=1:length(XSmt{1})
    % Matrix A coefficients
    temp = Param.X0 * Param.X0' + Param.W0;
    for  k=1:K-1
         temp = temp + XSmt{k}*XSmt{k}' + SSmt{k};
    end
    tempA = temp;
    % Matrix B coefficients
    if size(Uk,1)
        temp = zeros(size(Uk,2),size(Uk,2));
        for  k=1:K-1
            temp = temp + Uk(k,:)'*Uk(k,:);
        end
        tempB = temp;
    else
        tempB = [];
    end
    % Matrix AB coefficients
    if size(Uk,1)
        temp =  zeros(length(XSmt{1}),size(Uk,2));
        for  k=1:K-1
            temp = temp + XSmt{k}*Uk(k,:);
        end
        tempAB = temp;
    else
        tempAB = [];
    end
    % combine all three and take inverese
    temp = Param.Wk(r,r)* pinv([tempA   tempAB;
                               tempAB'  tempB]);
    COV_X{r}.A  = temp(1:size(Param.Ak,1),1:size(Param.Ak,1));
    COV_X{r}.AB = temp(1:size(Param.Ak),size(Param.Ak)+1:end);
    COV_X{r}.B  = temp(size(Param.Ak,1)+1:end,size(Param.Ak,1)+1:end);
    COV_X{r}.W  = Param.Wk(r,r)^2*2/K;
    COV_X{r}.SE_W = sqrt(diag(COV_X{r}.W));
    % we can also add covraince between A&W, B&W
end

%% Find Covariance for Continuous Part
% we assume fully observed data
xM = Param.xM;
% Normal
COV_C = [];
if DISTR(1)==1
    [MCk,MDk] = compass_Tk(In,Param);
    % covariance calculatuion
    tXSmt  = xM * XSmt{1};
    tSSmt  = xM * SSmt{1}* xM';
    tempA  = (MCk{1}'*MCk{1}).*(tXSmt*tXSmt'+tSSmt);
    tempB  = (MDk'*MDk).*(In(1,:)'*In(1,:));
    tempAB = (MCk{1}'.*tXSmt)*(MDk.*In(1,:));
    for  k=2:K
        tXSmt  = xM * XSmt{k};
        tSSmt  = xM * SSmt{k}* xM';
        
        tempA  = tempA  + (MCk{k}'*MCk{k}).*(tXSmt*tXSmt'+tSSmt);
        tempB  = tempB  + (MDk'*MDk).*(In(k,:)'*In(k,:));
        tempAB = tempAB + (MCk{k}'.*tXSmt)*(MDk.*In(k,:));
    end
    temp   = [tempA   tempAB;
            tempAB' tempB];
    ind    = [1:size(tempA,1) size(tempA,1)+find(MDk)];
    temp   = temp(:,ind); 
    temp   = temp(ind,:);
    temp_x = pinv(temp)*Param.Vk;
    temp   = zeros(length(MCk{1})+length(MDk));
    temp(ind,ind) = temp_x;
    % covariance   
    COV_C.C  = temp(1:length(Param.Ck),1:length(Param.Ck));
    COV_C.SE_C = sqrt(diag(COV_C.C));
    COV_C.CD = temp(1:length(Param.Ck),length(Param.Ck)+1:end);
    COV_C.D  = temp(length(Param.Ck)+1:end,length(Param.Ck)+1:end);
    COV_C.SE_D = sqrt(diag(COV_C.D));
    COV_C.V  = Param.Vk^2*2/K;
    COV_C.SE_V = sqrt(COV_C.V);
    % we can also add the C&V, D&V as well
end
% Gamma
if DISTR(1)== 2
    [MCk,MDk] = compass_Tk(In,Param);
    % model parameters
    Ck = Param.Ck;
    Dk = Param.Dk;
    Vk  = Param.Vk;
    S  = Param.S;
    % draw samples per trial
    N    = 1000;
    % map to larger space
    Xs   = xM * mvnrnd(XSmt{1},SSmt{1},N)';
    
    % calculate probability
    temp  = (Ck.*MCk{1})* Xs + (Dk.*MDk)* In(1,:)';
    Mx    = exp(-temp);
    
    Y     = Yn(1)-S;
    Xs_a  = repmat(Mx,length(Ck),1).*Xs;
    
    tempCC  = - Vk* Y *(MCk{1}'*MCk{1}).*(Xs_a * Xs')/N;
    tempCD  = - Vk* Y *(MCk{1}.*sum(Xs_a,2)'/N)'*(MDk.*In(1,:));
    tempDD  = - Vk* Y *(MDk'*MDk ).*(In(1,:)'*In(1,:))* sum(Mx)/N;
    
    tempCV  = - MCk{1}.*sum(Xs-Y*Xs_a,2)'/N;
    tempCS  = - Vk * MCk{1}.*sum(Xs_a,2)'/N;
    
    tempDV  = - MDk * (1-Y*sum(Mx)/N);
    tempDS  = - Vk* MDk * sum(Mx)/N;
    
    tempSS  = (1-Vk)/Y^2;
    tempVV  = - psi(1,Vk)+(1/Vk);
    tempSV  = (-1/Y) + sum(Mx)/N;
        
    for k = 2:K
        % map to larger space
        Xs   = xM * mvnrnd(XSmt{k},SSmt{k},N)';
        % calculate probability
        temp = (Ck.*MCk{k})* Xs + (Dk.*MDk)*In(k,:)';
        Mx   = exp(-temp);
        
        Y    = Yn(k)-S;
        Xs_a = repmat(Mx,length(Ck),1).*Xs;
       
        tempCC  = tempCC -Vk*Y*(MCk{k}'*MCk{k}).*(Xs_a * Xs')/N;
        tempDD  = tempDD-Vk*Y*(MDk'*MDk).*(In(k,:)'*In(k,:))*sum(Mx)/N;
        tempCD  = tempCD-Vk*Y*(MCk{k}.*sum(Xs_a,2)'/N)'*(MDk.*In(k,:));

        tempCV  = tempCV - MCk{k}.*sum(Xs-Y*Xs_a,2)'/N;
        tempCS  = tempCS - Vk*MCk{k}.*sum(Xs_a,2)'/N;

        tempDV  = tempDV -MDk  * (1-Y*sum(Mx)/N);
        tempDS  = tempDS -Vk* MDk * sum(Mx)/N;

        tempSS  = tempSS + (1-Vk)/Y^2;
        tempVV  = tempVV - psi(1,Vk)+(1/Vk);
        tempSV  = tempSV +(-1/Y) + sum(Mx)/N;

    end
    
    temp  = [tempCC    tempCD;
             tempCD'   tempDD];
    n_col = [tempCV tempDV]';
    temp  = [temp n_col];
    n_col = [tempCS tempDS]';
    temp  = [temp n_col];
    n_row = [tempCV tempDV tempVV tempSV];
    temp  = [temp;n_row];
    n_row = [tempCS tempDS tempSV tempSS];
    temp  = [temp;n_row];
    %% take out zero terms, calculate inverse, and then...
    if S~=0
        ind    = [1:size(tempCC,1) size(tempCC,1)+find(MDk)  size(tempCC,1)+length(MDk)+1 size(tempCC,1)+length(MDk)+2];
        temp   = temp(:,ind); 
        temp   = temp(ind,:);
        temp_x = -pinv(temp);
        temp   = zeros(length(MCk{1})+length(MDk)+2);
        temp(ind,ind)= temp_x;
    else
        ind    = [1:size(tempCC,1) size(tempCC,1)+find(MDk)  size(tempCC,1)+length(MDk)+1];
        temp   = temp(:,ind); 
        temp   = temp(ind,:);
        temp_x = -pinv(temp);
        temp   = zeros(length(MCk{1})+length(MDk)+2);
        temp(ind,ind)= temp_x;
    end
    % covariance   
    COV_C.C  = temp(1:length(Param.Ck),1:length(Param.Ck));
    COV_C.SE_C = sqrt(diag(COV_C.C));
    COV_C.D  = temp(length(Param.Ck)+1:length(Param.Ck)+length(Param.Dk),length(Param.Ck)+1:length(Param.Ck)+length(Param.Dk));
    COV_C.SE_D = sqrt(diag(COV_C.D));
    COV_C.S  = temp(end,end);
    COV_C.Vk  = temp(end-1,end-1);
    
    COV_C.CD = temp(1:length(Param.Ck),length(Param.Ck)+1:length(Param.Ck)+length(Param.Dk));
    COV_C.CV = temp(1:length(Param.Ck),length(Param.Ck)+length(Param.Dk)+1);
    COV_C.CS = temp(1:length(Param.Ck),length(Param.Ck)+length(Param.Dk)+2);
    COV_C.DV = temp(length(Param.Ck)+1:length(Param.Ck)+length(Param.Dk),length(Param.Ck)+length(Param.Dk)+1);
    COV_C.DS = temp(length(Param.Ck)+1:length(Param.Ck)+length(Param.Dk),length(Param.Ck)+length(Param.Dk)+2);
    COV_C.VS = temp(end,end-1);
    
end
%% Find Covariance for Discrete Part
% we assume fully observed data
COV_D = [];
if DISTR(2)==1
    [MEk,MFk] = compass_Qk(Ib,Param);
    % model parameters
    Ek   = Param.Ek;
    Fk   = Param.Fk;
    % draw samples per trial
    N    = 1000;
    % map to larger space
    Xs   = xM * mvnrnd(XSmt{1},SSmt{1},N)';
    % calculate probability
    temp = (Ek.*MEk{1})* Xs + (Fk.*MFk) * Ib(1,:)';
    Ps   = exp(temp)./(1+exp(temp));
    % duplicate X, one with P(1-P) as weight
    Xs_a = repmat(Ps.*(1-Ps),length(Ek),1).*Xs;
    Xs_b = Xs;
    % now calculate fisher
    tempA  = (MEk{1}'*MEk{1}).*(Xs_a*Xs_b')/N;
    tempB  = (MFk'*MFk).*(In(1,:)'*In(1,:))* sum(Ps.*(1-Ps))/N;
    tempAB = (MEk{1}.*sum(Xs_a,2)'/N)'*(MFk.*In(1,:));
    for k=2:K
        Xs   = xM * mvnrnd(XSmt{k},SSmt{k},N)';
        
        temp = (Ek.*MEk{k})*Xs+(Fk.*MFk)*Ib(k,:)';
        Ps   = exp(temp)./(1+exp(temp));
        Xs_a = repmat(Ps.*(1-Ps),length(Ek),1).*Xs;
        Xs_b = Xs;
    
        tempA  = tempA + (MEk{k}'*MEk{k}).*(Xs_a*Xs_b')/N;
        tempB  = tempB + (MFk'*MFk).*(In(k,:)'*In(k,:))*sum(Ps.*(1-Ps))/N;
        tempAB = tempAB+ (MEk{k}.*sum(Xs_a,2)'/N)'*(MFk.*In(k,:));
    end
    temp   = [tempA tempAB;
              tempAB' tempB];
    ind    = [1:size(tempA,1) size(tempA,1)+find(MFk)];
    temp   = temp(:,ind); 
    temp   = temp(ind,:);
    temp_x = pinv(temp);
    temp   = zeros(length(MEk{1})+length(MFk));
    temp(ind,ind)= temp_x;
    % covariance   
    COV_D.E  = temp(1:length(Param.Ek),1:length(Param.Ek));
    COV_D.SE_E = sqrt(diag(COV_D.E));
    COV_D.EF = temp(1:length(Param.Ek),length(Param.Ek)+1:end);
    COV_D.F  = temp(length(Param.Ek)+1:end,length(Param.Ek)+1:end);
    COV_D.SE_F = sqrt(diag(COV_D.F));
end
    

