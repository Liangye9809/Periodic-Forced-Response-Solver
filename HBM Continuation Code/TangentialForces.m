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



% function [T,w] = TangentialForce(xt,wt,kt,mu,FN)
%     Nx = size(xt,1); % xt is column vector
%     T = zeros(size(xt));
%     w = zeros(size(xt));
%     kt = kt(:); % if it is row vector, change to column vector
%     mu = mu(:);
%     FN = FN(:);
%     % if it is a 1*1 scale, change to Nx*1 column vector, same elements
%     if size(FN,1) == 1
%         FN = FN * ones(Nx,1);
%     end
%     if size(kt,1) == 1
%         kt = kt * ones(Nx,1);
%     end
%     if size(mu,1) == 1
%         mu = mu * ones(Nx,1);
%     end
% 
%     for i = 1:Nx
%         if FN(i) > 0
%             T(i) = kt(i) * (xt(i) - wt(i));
%             if abs(T(i)) < mu(i) * FN(i)
%                 w(i) = wt(i);
%             else
%                 sg = sign(T(i));
%                 T(i) = sg * mu(i) * FN(i);
%                 w(i) = xt(i) - sg * mu(i) * FN(i) / kt(i);
%             end
%         else
%             T(i) = 0;
%             w(i) = xt(i);
%         end
%     end
% 
% end

