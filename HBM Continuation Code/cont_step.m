function [x, lambda, tx, tlambda] = cont_step(func, jacob, derlambda, params)
%   此函数中，输入参数中包含x0,tx0,tlambda0，理论上需要先求解lambda0处的x0，然后再求解tx0和tlambda0。再通过ds求解下一个点的解
%   但是传递参数中包括了x0,tx0,tlambda0，此处x0是根据lambda0事先求解的，而tx0和tlambda0可以人为定义一个ds的方向，也可以通过求解(x0,lambda0)处的切线后输入。
x0 = params.cont.x0;
lambda0 = params.cont.lambda0;
ds = params.cont.ds;
tx0 = params.cont.tx0;
tlambda0 = params.cont.tlambda0;
epsx = params.Newton.epsx;
epsf = params.Newton.epsf;
maxiter = params.Newton.maxiter;
x = x0 + ds*tx0;
lambda = lambda0 + ds*tlambda0;
z = [x
     lambda];
% disp(':   k   x1   x2   lambda  errorx  errorf  ')
G = [func(x, lambda, params.func) 
     tx0'*(x-x0) + tlambda0*(lambda-lambda0)-ds];
for k=1:maxiter
    JG = [jacob(x, lambda, params.func) derlambda(x, lambda, params.func)
          tx0' tlambda0];
    dz = -JG\G;
    z = z + dz;
    x = z(1:end-1);
    lambda = z(end);
    G = [func(x, lambda, params.func) 
         tx0'*(x-x0) + tlambda0*(lambda-lambda0)-ds];
    errorx = norm(dz) / norm(z);
    errorf = norm(G);
    % format short g
    % disp([k z(1) z(2) lambda errorx  errorf])
    if (errorx < epsx) && (errorf < epsf)
        t0 = [tx0; tlambda0];
        t = null([jacob(x, lambda, params.func) derlambda(x, lambda, params.func)]);
        t = t * sign(t0'*t);
        tx = t(1:end-1);
        tlambda = t(end);
        format short g
        disp([params.cont.step  z(1) z(2) lambda  errorx  errorf   k])
        return
    end
end

end

