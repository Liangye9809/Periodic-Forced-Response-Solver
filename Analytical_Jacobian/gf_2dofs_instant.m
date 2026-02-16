function [F, w, Mf] = gf_2dofs_instant(x, kn, xn0, mu, kt, w_in, C_) 

   
    F = zeros(size(x));
    Mf = zeros(3, 1); % condition of friction matix
    w = w_in;
    
    
    F(2) = NormalForces(x(2), kn, xn0); % Fn
    [F(1), w, Mf(1:2)] = TangentialForces(x(1), w, kt, mu, F(2), C_); % T

    if abs(F(2)) > 0
        Mf(3) = 1;
    end

end

% C is condition of friction, [stick; slip] 1 is active, 0 is non.

% function [T, w, C] = TangentialForces(xt, wt, kt, mu, FN, C_)
%     if FN > 0
%         T = kt * (xt - wt);
%         if abs(T) < mu * FN
%             w = wt;
%             C = [1; 0]; % 100% stick
%         else
% 
%             if abs(T) > mu * FN 
%                 sg = sign(T);
%                 T = sg * mu * FN;
%                 w = xt - sg * mu * FN / kt;
%                 C = sg * [0; 1]; % 100% slip
%             else
%                 sg = sign(T);
%                 T = sg * mu * FN;
%                 w = xt - sg * mu * FN / kt;
%                 C = C_; % slip or stick, equal to previous time instant 
%             end
% 
%         end
%     else
%         T = zeros(size(xt));
%         w = xt;
%         C = [0; 0]; % gap, no stick nor slip
%     end
% end

% C is condition of friction, [stick; slip] 1 is active, 0 is non.

function [T, w, C] = TangentialForces(xt, wt, kt, mu, FN, C_)
    if FN > 0
        T = kt * (xt - wt);
        if abs(T) < mu * FN
            w = wt;
            C = [1; 0]; % 100% stick
        else
            sg = sign(T);
            T = sg * mu * FN;
            w = xt - sg * mu * FN / kt;

            C = sg * [0; 1]; 

            % C = C_; % slip or stick, equal to previous time instant 


        end
    else
        T = zeros(size(xt));
        w = xt;
        C = [0; 0]; % gap, no stick nor slip
    end
end


function FN = NormalForces(xn, kn, xn0)
    u = xn - xn0;
    FN = max(0, kn .* u);
end