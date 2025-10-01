function F = lp_mode_field(R, PHI, a, l, u, w)
% V jádře: J_l(u r / a) * cos(l*phi)
% V plášti: konst * K_l(w r / a), konst zajišťuje spojitost na r=a
    F = zeros(size(R));
    inside = R <= a; outside = ~inside;

    % jádro
    Jr = besselj(l, u * R(inside) / a);
    if l == 0
        ang_in = ones(size(Jr));
    else
        ang_in = cos(l * PHI(inside));
    end
    F(inside) = Jr .* ang_in;

    % plynulý přechod – škálování tak, aby F je spojité na r=a
    Fa = besselj(l, u);           % hodnota v r=a (J_l(u))
    Ka = besselk(l, w);           % hodnota v r=a (K_l(w))
    scale = Fa / Ka;              % aby K_l(w r/a) na r=a = J_l(u)

    Kr = besselk(l, w * R(outside) / a);
    if l == 0
        ang_out = ones(size(Kr));
    else
        ang_out = cos(l * PHI(outside));
    end
    F(outside) = scale * Kr .* ang_out;

    % normalizace (jen pro vzhled)
    F = F / max(abs(F(:)) + 1e-12);
end