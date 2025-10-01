% plot_lp_modes.m
% 9 LP modů ve vláknu (vizualizace intenzity / pole)
% Autor: ZFOT – MATLAB ukázka

clear; clc;

% ---------- Parametry vlákna a výpočtu ----------
a   = 1.0;          % poloměr jádra (libovolná normalizace)
V   = 8.0;          % normalizovaná frekvence (dostatečně velká pro dané mody)
N   = 401;          % rozlišení mřížky
xy  = 1.6*a;        % rozsah zobrazení [-xy, xy]

% Mřížka
x = linspace(-xy, xy, N);
y = linspace(-xy, xy, N);
[XX, YY] = meshgrid(x, y);
RR  = hypot(XX,YY);
PHI = atan2(YY,XX);

% Seznam modů [l, m] v pořadí pro 3x3 panel
modes = [ 0 1; 1 1; 2 1; 3 1;
          0 2; 1 2; 2 2; 3 2;
          0 3; 1 3; 2 3; 3 3;
          0 4; 1 4; 2 4; 3 4];

% ---------- Vykreslení ----------
figure('Color','w','Position',[100 100 980 880]);

for k = 1:size(modes,1)
    l = modes(k,1); m = modes(k,2);

    % kořen Besselovy funkce J_l pro dané (l,m) – aproximační start a fzero
    alpha = bessel_zero(l, m);

    % rozdělení V: u = alpha, w = sqrt(V^2 - u^2) (pro vizualizaci volíme V>alpha)
    u = alpha;
    w = sqrt(max(V^2 - u^2, 1e-6));  % ochrana proti zápornému kvůli zaokrouhlení

    % pole (neškálované) – jedna z ortogonálních orientací (cos l*phi)
    F = lp_mode_field(RR, PHI, a, l, u, w);

    % zobrazení intenzity |F|^2 nebo real(F) dle preferencí:
    I = abs(F).^2;  % intenzita
    %I = real(F);   % alternativně pole s kladnou/zápornou fází

    subplot(4,4,k);
    imagesc(x/a, y/a, I); axis image off;


    title(sprintf('LP_{%d%d}', l, m), 'FontWeight','bold');
    set(gca,'YDir','normal');
    % colormap(gca, blueWhiteRed(1024)); % příjemná divergující mapa
    colormap(gca, hot);
    % Normalizace kontrastu pro čitelnost
    clim = prctile(I(:), [1 99.5]);
    set(gca,'CLim',clim);
end

t = sgtitle('LP módy (normalizované souřadnice r/a)','FontWeight','bold'); %#ok<NASGU>


