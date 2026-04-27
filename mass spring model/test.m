load('data\Analytical Jacobian results 2.0 simple sin\para_gap_to_stick_point.mat');
n = size(para.x_cont, 2);

for i = 1:n
    xt = Fourier_to_Time(para.x_cont(:, i), )
end



function xct = Fourier_to_Time(Xc, H, Nx, E)

    for i = 1:3 * Nx
        r1 = (2 * H + 1) * (i - 1) + 1;
        r2 = (2 * H + 1) * i;
        X(:,i) = Xc(r1:r2); % reorder in dofs in column
    end
    xct = E * X; 

end