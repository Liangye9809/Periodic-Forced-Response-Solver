function F = perform_aft(Xc, H, E, EH, fc)
    
    X1 = Xc(1:(2 * H + 1));
    X2 = Xc((2 * H + 2):(4 * H + 2));
    X3 = Xc((4 * H + 3):end);

   
    x1 = E * X1;
    x2 = E * X2;
    x3 = E * X3;

    N = size(x1, 1);
    N_cycles = 2;
    nlforce = zeros(N, 3);
    
    for j = 1:N_cycles
        for i = 1:N
            ggStruct = gg(x1(i), x2(i), x3(i), fc);
            nlforce(i, :) = ggStruct.F;
            fc.w = ggStruct.w;
        end
    end

    F.F = EH * nlforce;
    F.w = ggStruct.w;
end


    