function [F, w, Mf, dxdn] = gf_2dofs_instant(x, kn, xn0, mu, kt, w_in, C_, xn) 

   
    F = zeros(size(x));
    Mf = zeros(3, 1); % condition of friction matix
    w = w_in;
    
    
    F(2) = NormalForces(x(2), kn, xn0); % Fn
    [F(1), w, Mf(1:2), dxdn] = TangentialForces(x(1), w, kt, mu, F(2), C_, xn); % T

    if abs(F(2)) > 0
        Mf(3) = 1;
    end

end


% C is condition of friction, [stick; slip] 1 is active, 0 is non.
% corrected Tangential force, interpolation w

% function [T, w, C, dxdn] = TangentialForces(xt, wt, kt, mu, FN, C_, xn)
%     dxdn = [0, 0];
%     wp = wt;
%     if FN > 0
%         T = kt * (xt - wp);
%         if xn(1) <= 0 % previous is gap, means gap to contact
%             dxdti = (xt - wt) / (xn(2) - xn(1) + 1e-16);
%             wp = wt - xn(1) * dxdti;
%             T = kt * (xt - wp); 
%             % dxdn = [1, dxdti];
%         end
%         if abs(T) < mu * FN
%             w = wp; 
%             C = [1; 0]; % stick
%             if xn(1) <= 0 % gap to stick
%                 dxdn = [1, dxdti]; % gap to stick correction part
%             end
%         else
%             sg = sign(T);
%             T = sg * mu * FN;
%             w = xt - sg * mu * FN / kt;
%             C = sg * [0; 1]; % 100% slip
%         end
%     else
%         T = zeros(size(xt));
%         w = xt;
%         C = [0; 0]; % gap, no stick nor slip
%     end
% end


% C is condition of friction, [stick; slip] 1 is active, 0 is non.
% no interpolation of w

function [T, w, C, dxdn] = TangentialForces(xt, wt, kt, mu, FN, C_, xn)
    dxdn = [0, 0];
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN % stick
            w = wt; 
            C = [1; 0]; % stick
            if xn(1) <= 0 % previous is gap, means gap to stick
                dxdti = (xt - wt) / (xn(2) - xn(1) + 1e-16);
                dxdn = [1, dxdti]; % get this value for analytical Jacobian
            end
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


% function [T, w, C, dxdn] = TangentialForces(xt, wt, kt, mu, FN, C_, xn)
%     dxdn = [0, 0];
%     if FN > 0
%         T = kt * (xt - wt);
%         if abs(T) <= mu * FN
%             w = wt; % xt⁻
%             if and(C_ == [0; 0], abs(1 - xn(2) / xn(1)) > 1e-5)% previous is gap, means gap to stick
%                 dxdti = (xt - wt) / (xn(2) - xn(1));
%                 w = wt - xn(1) * dxdti;
%                 T = kt * (xt - w); % Does it need to be compare with mu*Fn?
%                 dxdn = [1, dxdti];
%             end
%             C = [1; 0]; % stick
%         else
%             sg = sign(T);
%             T = sg * mu * FN;
%             w = xt - sg * mu * FN / kt;
%             C = sg * [0; 1]; % 100% slip
%         end
%     else
%         T = zeros(size(xt));
%         w = xt;
%         C = [0; 0]; % gap, no stick nor slip
%     end
% end