
function [T, w] = TangentialForces(xt, wt, kt, mu, FN)
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN
            w = wt;
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;
        end
    else
        T = zeros(size(xt));
        w = xt;
    end
end

