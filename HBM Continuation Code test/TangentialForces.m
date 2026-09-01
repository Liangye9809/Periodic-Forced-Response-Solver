% function [T, w, Flag, dxdn] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre)
%     dxdn = [0, 1];
%     if FN > 0
%         T = kt * (xt - wt);
%         if abs(T) < mu * FN
%             w = wt;
%             Flag = 2; % stick
%             if xn_pre <= 0 % previous is gap, means gap to stick
%                 dxdti = (xt - wt) / (xn - xn_pre);
%                 w = wt - xn_pre * dxdti; % update w
%                 T = kt * (xt - w); % Does it need to be compare with mu*Fn?
%                 dxdn = [1, dxdti];
%             end
%         else
%             sg = sign(T);
%             T = sg * mu * FN;
%             w = xt - sg * mu * FN / kt;
%             Flag = sg; % slip
%         end
%     else
%         T = zeros(size(xt));
%         w = xt;
%         Flag = 0; % gap
%     end
% end



function [T, w, Flag, dxdn] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre, kn, dxt, dxn, Flag_pre, dtheta)
    dxdn = [0, 1];
    wp = wt; % previous w
    if FN > 0 % contact
        T = kt * (xt - wp);
        if abs(Flag_pre) == 1
            q = dxt - Flag_pre * mu * kn / kt * dxn; % dw/dtheta on the slip branch
            Q = Flag_pre * q;
            if Q(1) > 0 && Q(2) <= 0 && (abs(T) < mu * FN || sign(T) == Flag_pre)
                alpha = Q(1) / (Q(1) - Q(2));
                wp = wt + 0.5 * dtheta * alpha * q(1); % integral of linear q up to q(theta*) = 0
                T = kt * (xt - wp);
                w = wp;
                Flag = 2;
                return
            end
        end
        if xn_pre <= 0 % previous is gap, means gap to stick
            dxdti = (xt - wt) / ((xn - xn_pre) + 1e-16);
            wp = wt - xn_pre * dxdti; % update w
            T = kt * (xt - wp);
            dxdn = [1, dxdti];
        end
        if abs(T) < mu * FN
            w = wp;
            Flag = 2; % stick
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

% function [T, w, Flag, dxdn] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre)
%     if FN > 0
%         T = kt * (xt - wt);
%         if abs(T) < mu * FN
%             w = wt; 
%             Flag = 2;
%         else
%             sg = sign(T);
%             T = sg * mu * FN;
%             w = xt - sg * mu * FN / kt;
%             Flag = sg;
%         end
%     else
%         T = zeros(size(xt));
%         w = xt;
%         Flag = 0;
%     end
%     dxdn = 0;
% end
