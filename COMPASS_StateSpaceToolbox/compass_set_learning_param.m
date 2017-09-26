function Param = compass_set_learning_param(Param,Iter,UpdateStateParam,UpdateStateNoise,UpdateStateX0,UpdateCModelParam,UpdateCModelNoise,UpdateDModelParam,DiagonalA,UpdateMode,UpdateCModelShift)

    %% Set Learning Rule
    % number of iteration
    Param.Iter = Iter;
    % Hidden model update
    % we can set to 0, which means neither A or B is being updated
    % we can set to 1, which means both A and B are being updated
    % we can set to 2, which means A is being fixed and B is getting updated
    Param.UpdateStateParam =  UpdateStateParam;
    
    Param.UpdateStateNoise =  UpdateStateNoise;
    Param.UpdateStateX0    =  UpdateStateX0;
    % continuous model update
    Param.UpdateCModelParam = UpdateCModelParam;
    Param.UpdateCModelNoise = UpdateCModelNoise;
    % discrete model update
    Param.UpdateDModelParam = UpdateDModelParam;
    % marix A
    Param.DiagonalA = DiagonalA;
    % Check Update Model
    Param.UpdateMode=UpdateMode;
    % Shift in Gamma
    Param.UpdateCModelShift=UpdateCModelShift;

end

