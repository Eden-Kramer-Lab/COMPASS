function Param = compass_create_state_space(nx,nUk,nIn,nIb,xM,cLink,cLinkUpdate,dLink,dLinkUpdate)
    
    %% Build The State-Space Model
    nc = 1;
    nd = 1;
    %% Dimension of Components
    % Size of X vector
    Param.nx = nx;   
    Param.nc = 1;
    Param.nd = 1;
    
    %% Input and Their Link Functions
    % Input and its link to the continuous model
    Param.nIn = nIn;         % Bias, I, V, I2C, C2I
    
    if ~isempty(cLink) && ~isempty(cLinkUpdate)
        % The cLinkMap - continuous link map
        Param.cLinkMap = zeros(nc,size(xM,1));
        for i=1:nc
            Param.cLinkMap(i,:)   = cLink(i,:);    
            Param.cLinkUpdate(i,:)= cLinkUpdate(i,:);    
        end
    end
    % Input and its link to the discrete model
    Param.nIb = nIb;         % Bias, I, V, I2C, C2I
    if ~isempty(dLink) && ~isempty(dLinkUpdate)
        % The dLinkMap - discrete link map
        Param.dLinkMap = zeros(nd,size(xM,1));
        for i=1:nd
            Param.dLinkMap(i,:)   = dLink(i,:); 
            Param.dLinkUpdate(i,:)= dLinkUpdate(i,:);    
        end
    end

    %% Model Paramaters
    % Hidden state model - 
    Param.Ak = eye(Param.nx,Param.nx)*1;
    Param.Bk = ones(Param.nx,nUk)*0.0;      % input
    Param.Wk = eye(Param.nx,Param.nx)*0.05; % iid noise
    Param.X0 = zeros(Param.nx,1);           % initial x0 is set to 0
    Param.W0 = eye(Param.nx,Param.nx)*1;  % initial x0 is set to 0

    % Continuous model
    Param.Ck = ones(nc,size(xM,1));      % coefficients of the x
    Param.Dk = ones(nc,Param.nIn);     % input parameters
    % we need to drop some input from update
    Param.cConstantUpdate= ones(nc,Param.nIn);    
    if ~isempty(cLink) && ~isempty(cLinkUpdate)
        for i=1:nc
            ind = find(Param.cLinkMap(i,:));
            Param.cConstantUpdate(i,Param.cLinkMap(i,ind))=0;
        end
    end
    Param.Vk = eye(nc,nc)*0.01; % noise


    % Disceret model
    Param.Ek = ones(nd,size(xM,1));      % coefficients of the x
    Param.Fk = ones(nd,Param.nIb);     % input parameters
    % we need to drop some input from update
    if ~isempty(dLink) && ~isempty(dLinkUpdate)
        Param.dConstantUpdate= ones(nd,Param.nIb);   
        for i=1:nd
            ind = find(Param.dLinkMap(i,:));
            Param.dConstantUpdate(i,Param.dLinkMap(i,ind))=0;
        end
    end
   
    % xM is the X Mapping
    Param.xM = xM;
    
    % two extra paramaters for Gamma (Param.Vk is beging treated as Dispression)
    Param.S  = 0;
    
    % set censor_time + processing model
    Param.censor_time = 1;
    Param.censor_mode = 1;  % mode 1: sampling, mode 2: full likelihood
    % mode 1 is applicable on all data types
    % mode 2 is defined on continuous variable
    Param.censor_update_mode=1; % how the run the filter rule
    
    % it define how the Kalman Filter is getting updated - two different
    % mehtods
    Param.UpdateMode=1;

end

