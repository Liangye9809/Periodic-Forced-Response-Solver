function F = g(xt) % xt, each rows correspond each time, columns are different dofs
%% sqare
% F = xt .^ 2;

%%
x1 = xt(:, 1);
x2 = xt(:, 2);
x3 = xt(:, 3);

F1 = x1 .* x2 .* x3;
F2 = exp(x1 + x3) .* x2;
F3 = 1 ./ (1 + x1.^2 + x2.^2 + x3.^2); 

F = [F1, F2, F3];


end