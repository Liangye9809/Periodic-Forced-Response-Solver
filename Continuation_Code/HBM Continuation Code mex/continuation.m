function [x_cont, omega_cont] = continuation(FUNC, JACOB, JLAMBDA, params)
    
    % First, Calculate the initial point, with ds = 0, 
    ds = params.cont.ds;
    omega_0 = params.cont.omega_0;
    omega_end = params.cont.omega_end;
    % calculate the initial point
    params.cont.ds = 0;
    params.cont.step = 1;
    disp('continuation information');
    disp('step:   x1   x2   lambda  errorx  errorf     kmax');
    [x,omega,tx,tomega,w] = cont_step(FUNC, JACOB, JLAMBDA, params);
    params.func.fc.w = w; % update w
    x_cont = x;
    omega_cont = omega;
    if omega_0 == omega_end
        return;
    end
    % calculate the rest points
    params.cont.ds = ds;
    params.cont.x0 = x;
    params.cont.tx0 = tx;
    params.cont.omega_0 = omega;
    params.cont.tomega0 = tomega;
    maxstep = params.cont.maxstep;
    for n = 2:maxstep
        params.cont.step = n;
        [x,omega,tx,tomega,w] = cont_step(FUNC, JACOB, JLAMBDA, params);
        params.func.fc.w = w; % update w
        params.cont.x0 = x;
        params.cont.omega_0 = omega;
        params.cont.tx0 = tx;
        params.cont.tomega0 = tomega;
        x_cont = [x_cont, x];
        omega_cont = [omega_cont, omega];
        if sign(omega_end - omega_0) * (omega + ds*tomega - omega_end) > 0 % the final point
            params.cont.ds = 0;
            params.cont.tx0 = zeros(size(tx));
            params.cont.tomega0 = sign(omega_end - omega_0) * 1;
            params.cont.omega_0 = omega_end;
            params.cont.step = n + 1;
            [x,omega,tx,tomega,w] = cont_step(FUNC, JACOB, JLAMBDA, params);
            params.func.fc.w = w; % update w
            x_cont = [x_cont, x];
            omega_cont = [omega_cont, omega];
            break;
        end
        
    end
end