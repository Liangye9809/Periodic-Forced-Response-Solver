% function [F, w] = g(xt, fc) % xt, each rows correspond each time, columns are different dofs
% %% cubic 
% % F = 1 .* 1e8 .* xt.^3;
% % w = fc.w;
% 
% %% friction
% [F, w] = g_mex(xt, fc);
% 
% end


function [Ft, wt, flag] = g(xt, kn, xn0, mu, kt, w_in, nloop) % xt, each rows correspond each time, columns are different dofs

    [N, M] = size(xt);
    Nx = M / 3;

    Ft = zeros(nloop * N, M);
    wt = zeros(2, Nx, nloop * N);

    flag = zeros(2, Nx, nloop * N);
    xnt = xt(:, 3:3:end);
    

    for j = 1:nloop
        for i = 1:N
            ipre = mod(i - 2, N) + 1;
            [Ft((j - 1) * N + i, :), wtemp, flagi] = gf(xt(i, :), kn, xn0, mu, kt, w_in, xnt(ipre, :));
            w_in = wtemp; 
            wt(:, :, (j - 1) * N + i) = wtemp;
            flag(:, :, (j - 1) * N + i) = flagi;
            
        end
    end

    
end
