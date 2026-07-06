function [x_cont, omega_cont, k_cont, w_cont, stick_cont, slipP_cont, slipM_cont, gap_cont] = continuation(FUNC, JACOB, JLAMBDA, params)
    
    % First, Calculate the initial point, with ds = 0, 
    ds = params.cont.ds;
    omega_0 = params.cont.omega_0;
    omega_end = params.cont.omega_end;
    % calculate the initial point
    params.cont.ds = 0;
    params.cont.step = 1;
    disp('continuation information');
    disp('step:   x1   x2   lambda  errorx  errorf     kmax');
    [x, omega, tx, tomega, w, k, FlagState, stop] = cont_step(FUNC, JACOB, JLAMBDA, params);
    params.func.fc.w = w; % update w
    x_cont = x;
    k_cont = k;
    w_cont = w;
    stick_cont = FlagState(:, 1);
    slipP_cont = FlagState(:, 2);
    slipM_cont = FlagState(:, 3);
      gap_cont = FlagState(:, 4);
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
        [x, omega, tx, tomega, w, k, FlagState, stop] = cont_step(FUNC, JACOB, JLAMBDA, params);
        if stop == 1
            warning('iteration reaches the maximum step, continuation not converge');
            break;
        end
        params.func.fc.w = w; % update w
        params.cont.x0 = x;
        params.cont.omega_0 = omega;
        params.cont.tx0 = tx;
        params.cont.tomega0 = tomega;
        x_cont = [x_cont, x];
        omega_cont = [omega_cont, omega];
        k_cont = [k_cont, k];
        w_cont = [w_cont, w];
        stick_cont = [stick_cont, FlagState(:, 1)];
        slipP_cont = [slipP_cont, FlagState(:, 2)];
        slipM_cont = [slipM_cont, FlagState(:, 3)];
          gap_cont = [  gap_cont, FlagState(:, 4)];
        if sign(omega_end - omega_0) * (omega + ds*tomega - omega_end) > 0 % the final point
            params.cont.ds = 0;
            params.cont.tx0 = zeros(size(tx));
            params.cont.tomega0 = sign(omega_end - omega_0) * 1;
            params.cont.omega_0 = omega_end;
            params.cont.step = n + 1;
            [x, omega, tx, tomega, w, k, FlagState, stop] = cont_step(FUNC, JACOB, JLAMBDA, params);
            if stop == 1
                warning('iteration reaches the maximum step, continuation not converge');
                break;
            end
            params.func.fc.w = w; % update w
            x_cont = [x_cont, x];
            omega_cont = [omega_cont, omega];
            k_cont = [k_cont, k];
            w_cont = [w_cont, w];
            stick_cont = [stick_cont, FlagState(:, 1)];
            slipP_cont = [slipP_cont, FlagState(:, 2)];
            slipM_cont = [slipM_cont, FlagState(:, 3)];
              gap_cont = [  gap_cont, FlagState(:, 4)];
            break;
        end
        
    end
end