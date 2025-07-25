function F = force_fourier(X, xp, gxp, H, E, EH, fc)
    
    Xc = X;
    Xc(1) = Xc(1) + 2*xp(1);
    Xc(2*H + 2) = Xc(2*H + 2) + 2*xp(2);
    Xc(4*H + 3) = Xc(4*H + 3) + 2*xp(3);

    GcStruct = perform_aft(Xc, H, E, EH, fc);

    Gc = GcStruct.F;
    Gc(1) = Gc(1) - 2*gxp(1);
    Gc(2*H + 2) = Gc(2*H + 2) - 2*gxp(2);
    Gc(4*H + 3) = Gc(4*H + 3) - 2*gxp(3);

    F.F = Gc;
    F.w = GcStruct.w;
end