function [x, omega, tx, tomega, w] = cont_step(func, jacob, deromega, params)
x0 = params.cont.x0;
omega_0 = params.cont.omega_0;
ds = params.cont.ds;
tx0 = params.cont.tx0;
tomega0 = params.cont.tomega0;
epsx = params.Newton.epsx;
epsf = params.Newton.epsf;
maxiter = params.Newton.maxiter;
x = x0 + ds*tx0;
omega = omega_0 + ds*tomega0;
z = [x
     omega];
% disp(':   k   x1   x2   lambda  errorx  errorf  ')
% FUNCstruct = func(x, omega, params.func);
[F, w] = func(x, omega, params.func);
% params.func.fc.w = FUNCstruct.w; % update w
params.func.fc.w = w; % update w
G = [F 
     tx0'*(x-x0) + tomega0*(omega-omega_0)-ds];
for k=1:maxiter
    JG = [jacob(x, omega, params.func) deromega(x, omega, params.func)
          tx0' tomega0];
    A = [jacob(x, omega, params.func) deromega(x, omega, params.func)];
    e_ = svd(A);
    dz = -JG\G;
    z = z + dz;
    x = z(1:end-1);
    omega = z(end);
    % FUNCstruct = func(x, omega, params.func);
    [F, w] = func(x, omega, params.func);
    % params.func.fc.w = FUNCstruct.w; % update w
    params.func.fc.w = w; % update w
    G = [F 
         tx0'*(x-x0) + tomega0*(omega-omega_0)-ds];
    errorx = norm(dz) / norm(z);
    errorf = norm(G);
    % format short g
    % disp([k z(1) z(2) lambda errorx  errorf])
    if (errorx < epsx) && (errorf < epsf)
        t0 = [tx0; tomega0];
        t = null([jacob(x, omega, params.func) deromega(x, omega, params.func)]);
        t = t * sign(t0'*t);
        tx = t(1:end-1);
        tomega = t(end);
        format short g
        disp([params.cont.step  z(1) z(2) omega  errorx  errorf   k])
        % w = FUNCstruct.w; % update w
        return
    end
end


end

