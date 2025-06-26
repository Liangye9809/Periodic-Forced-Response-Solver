
omega02 = CB.CBmods.Kaa(1);
alpha = abs(CB.CB_F.Fa(1)) / omega02;
beta = abs(mu(1) * max(abs(Rx)) / kt(1));

CB.CB_F.Fa = CB.CB_F.Fa / (alpha * omega02);
CB.CB_F.Fx = CB.CB_F.Fx * beta / (alpha^2 * omega02);

Rx = Rx * beta / (alpha^2 * omega02);

CB.CB_MK.Max = CB.CB_MK.Max * (beta / alpha);
CB.CB_MK.Mxx = CB.CB_MK.Mxx * (beta / alpha)^2;

CB.CB_MK.Kaa = CB.CB_MK.Kaa / omega02;
CB.CBmods.Kaa = CB.CB_MK.Kaa;
CB.CB_MK.Kxx = CB.CB_MK.Kxx * (beta / alpha)^2 / omega02;

kn = kn * (beta / alpha)^2 / omega02;
kt = kt * (beta / alpha)^2 / omega02;

xi = xi * sqrt(omega02);

params.func.Nondimention.beta = beta;
params.func.Nondimention.alpha = alpha;
