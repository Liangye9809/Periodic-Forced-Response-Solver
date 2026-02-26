function [T, w, C] = TangentialForces(xt, wt, kt, mu, FN)
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) <= mu * FN
            w = wt;
            C = [1; 0]; % stick
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;
            C = sg * [0; 1]; % 100% slip
        end
    else
        T = zeros(size(xt));
        w = xt;
        C = [0; 0]; % gap, no stick nor slip
    end
end