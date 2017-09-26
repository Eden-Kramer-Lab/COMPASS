function Param = compass_set_censor_threshold_proc_mode(Param,censor_thr,censor_mode,update_mode)
    
    Param.censor_time = censor_thr;
    Param.censor_mode = censor_mode;  % mode 1: sampling, mode 2: full likelihood
    Param.censor_update_mode = update_mode; % eithe 1 or 2
    % mode 2 is valid when there is a continuous variable

end

