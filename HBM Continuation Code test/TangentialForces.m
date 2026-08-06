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



%% original logic
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


%% fix w in contact to gap
% function [T, w, Flag, dxdn] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre)
%     dxdn = [0, 1];
%     wp = wt; % previous w
%     if FN > 0 % contact
%         T = kt * (xt - wp);
%         if xn_pre <= 0 % previous is gap, means gap to stick
%             dxdti = (xt - wt) / ((xn - xn_pre) + 1e-16);
%             wp = wt - xn_pre * dxdti; % update w
%             T = kt * (xt - wp);
%             dxdn = [1, dxdti];
%         end
%         if abs(T) < mu * FN
%             w = wp;
%             Flag = 2; % stick
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

%% fix w contact-gap and slip-stick
function [T, w, Flag, dxdn] = TangentialForces(xt, wt, kt, mu, FN, xn, xn_pre, dxtdt_m, dxtdt_p, dxndt_m, dxndt_p, flag_m, N, kn)
    dxdn = [0, 1];
    wp = wt; % previous w
    if FN > 0 % contact
        T = kt * (xt - wp);
        if xn_pre <= 0 % previous is gap, means gap to stick
            dxdti = (xt - wt) / ((xn - xn_pre) + 1e-16);
            wp = wt - xn_pre * dxdti; % update w
            T = kt * (xt - wp);
            dxdn = [1, dxdti];
        end
        if abs(T) < mu * FN % stick
            if ismember(flag_m, [1,-1]) % slip to stick
                    dwdt_m = dxtdt_m - flag_m * mu * kn / kt * dxndt_m;
                    dwdt_p = dxtdt_p - flag_m * mu * kn / kt * dxndt_p;
                if dwdt_p * dwdt_m < 0 % different signs
                    fdt = dwdt_m / (dwdt_m - dwdt_p);
                    dt = 2 * pi / N;
                    % k_w = dxtdt_m + 0.5 * fdt * (dxtdt_p - dxtdt_m)...
                    %       - flag_m * mu * kn / kt * (dxndt_m + 0.5 * fdt * (dxndt_p - dxndt_m));
                    k_w = 0.5 * dwdt_m;
                    w = wp +  k_w * dt * fdt;
                    T = kt * (xt - w);
                else
                    w = wp;
                end
            else
                w = wp;
            end
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