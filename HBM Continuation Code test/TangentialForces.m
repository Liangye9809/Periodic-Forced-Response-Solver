function [T, w, Flag] = TangentialForces(xt, wt, kt, mu, FN, xn)
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN
            w = wt;
            Flag = 2; % stick
            if xn(1) <= 0 % previous is gap, means gap to stick
                dxdti = (xt - wt) / (xn(2) - xn(1));
                w = wt - xn(1) * dxdti; % update w
                T = kt * (xt - w); % Does it need to be compare with mu*Fn?
                dxdn = [1, dxdti];
            end
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;
            Flag = sg; % slip
        end
    else
        T = zeros(size(xt));
        w = xt;
        Flag = 0; % gap
    end
end