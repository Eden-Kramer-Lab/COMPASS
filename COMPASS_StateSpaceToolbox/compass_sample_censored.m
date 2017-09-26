function Yn = compass_sample_censored(DISTR,In,Ct,Dt,Param,Xm,Sm)
    
    if DISTR == 2
        % calculate mean at Xm
        Mx = exp(Ct * Xm + Dt * In');
        % calculate q2, exp(1/M)
        q2 = exp(-Dt*In')* exp(- Ct*Xm + 0.5 * Ct * Sm * Ct');

        % L0, L1, L2
        a  = Param.Vk;
        b  = Mx/Param.Vk;
        L0 = gamcdf(Param.censor_time-Param.S,a,b,'upper')*b^a*gamma(a);
        L1 = gamcdf(Param.censor_time-Param.S,a+1,b)*b^(a+1)*gamma(a+1);
        L2 = gamcdf(Param.censor_time-Param.S,a+2,b)*b^(a+2)*gamma(a+2);

        q3 = log(L0)+0.5*((L0*L2-L1*L1)/L0^2)*(Ct*Sm*Ct');


        %% call optimization function
        options = optimoptions('lsqnonlin','Display','off','DiffMaxChange',100,'MaxIter',1000);
        lower_bound = Param.censor_time;
        upper_bound = Inf;
        % model terms
        q1 = Param.Vk;
        % q2;
        % q3;
        % call optimzation function
        p0   = Param.censor_time+1000*eps;   
        Yn   = lsqnonlin(@OptimPoint,p0,lower_bound,upper_bound,options);
        % optimzization funcion
       
    else
        ok=1;
    end
    
     function f = OptimPoint(p)
            f  = (-(q1-1)*log(p)+q1*q2*p+q3)^2;
     end
    
end