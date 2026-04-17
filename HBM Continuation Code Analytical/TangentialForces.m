function [T, w, Flag, dxdn] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre)
    dxdn = [0, 1];
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN
            w = wt;
            Flag = 2; % stick
            if xn_pre <= 0 % previous is gap, means gap to stick
                dxdti = (xt - wt) / (xn - xn_pre);
                w = wt - xn_pre * dxdti; % update w
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