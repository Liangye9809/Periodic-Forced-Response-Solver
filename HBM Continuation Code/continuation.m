function [x_cont, lambda_cont] = continuation(FUNC, JACOB, JLAMBDA, params)
    % 先计算初值点的解，不调用cont_step函数，或者传入的ds = 0
    % 然后存储第一个点，按照正常的逻辑调用cont_step函数计算之后的点
    

    ds = params.cont.ds;
    lambda0 = params.cont.lambda0;
    lambda_end = params.cont.lambda_end;
    % calculate the initial point
    params.cont.ds = 0;
    params.cont.step = 1;
    disp('continuation information');
    disp('step:   x1   x2   lambda  errorx  errorf     kmax');
    [x,lambda,tx,tlambda] = cont_step(FUNC, JACOB, JLAMBDA, params);
    x_cont = x;
    lambda_cont = lambda;
    if lambda0 == lambda_end
        return;
    end
    % calculate the rest points
    params.cont.ds = ds;
    params.cont.x0 = x;
    params.cont.tx0 = tx;
    params.cont.lambda0 = lambda;
    params.cont.tlambda0 = tlambda;
    maxstep = params.cont.maxstep;
    for n = 2:maxstep
        params.cont.step = n;
        [x,lambda,tx,tlambda] = cont_step(FUNC, JACOB, JLAMBDA, params);
        params.cont.x0 = x;
        params.cont.lambda0 = lambda;
        params.cont.tx0 = tx;
        params.cont.tlambda0 = tlambda;
        x_cont = [x_cont, x];
        lambda_cont = [lambda_cont, lambda];
        if sign(lambda_end - lambda0) * (lambda + ds*tlambda - lambda_end) > 0 % the final point
            params.cont.ds = 0;
            params.cont.tx0 = zeros(size(tx));
            params.cont.tlambda0 = sign(lambda_end - lambda0) * 1;
            params.cont.lambda0 = lambda_end;
            params.cont.step = n + 1;
            [x,lambda,tx,tlambda] = cont_step(FUNC, JACOB, JLAMBDA, params);
            x_cont = [x_cont, x];
            lambda_cont = [lambda_cont, lambda];
            break;
        end
        % if omega > 1.05
        %     stophere = 1;
        % end
    end
end