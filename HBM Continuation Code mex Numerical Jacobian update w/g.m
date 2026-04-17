function [F, w] = g(xt, fc) % xt, each rows correspond each time, columns are different dofs
%% cubic 
% F = 1 .* 1e8 .* xt.^3;
% w = fc.w;

%% friction
[F, w] = g_mex(xt, fc);
    
end