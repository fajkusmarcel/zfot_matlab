%% Šíření EM vlny přes rozhraní dvou prostředí (kolmý dopad)
% y(z,t) ~ E_t(z,t), 1D, rovinná vlna
% - Zleva: prostředí n1 (incidentní + volitelně odražená vlna)
% - Zprava: prostředí n2 (přenesená vlna)
% - f je stanovena ve vakuu (λ0 = c/f), v prostředích platí v = c/n, λ = λ0/n
%
% Předpoklady: μ_r = 1 (nemagnetická média), Fresnelovy vztahy pro E při kolmém dopadu:
%   r = (n1 - n2) / (n1 + n2)              (amplitudový odraz E)
%   t = 2 n1 / (n1 + n2)                   (amplitudový přenos E)
%   R = |r|^2,  T = (n2/n1) |t|^2          (energetické poměry; R + T = 1)

clear; clc; close all;

%% --- Parametry simulace (můžeš volně měnit) --------------------------------
c   = 1;          % rychlost světla ve vakuu (normalizovaně c=1)
f   = 3.00;       % frekvence ve vakuu [Hz] (společná ve všech prostředích)
A0  = 1.00;       % amplituda incidentní vlny (E) v n1
phi = 0.00;       % počáteční fáze [rad]

n1  = 2.00;       % index lomu prostředí 1 (z<0)
n2  = 1.00;       % index lomu prostředí 2 (z>=0)

showReflection = false;   % =true: zobrazí se i odražená vlna v n1; =false: jen incident + transmit

% Časování
n_steps = 240;    % počet kroků animace
dt      = 0.01;   % časový krok [s]
Tpause  = 40;     % pauza mezi snímky [ms]

%% --- Odvozené veličiny ------------------------------------------------------
lambda0 = c / f;          % vlnová délka ve vakuu
omega   = 2*pi*f;         % úhlová frekvence
k0      = 2*pi / lambda0; % vlnové číslo ve vakuu

v1 = c / n1;    lambda1 = lambda0 / n1;   k1 = n1 * k0;
v2 = c / n2;    lambda2 = lambda0 / n2;   k2 = n2 * k0;

% Fresnelovy koeficienty (amplitudové pro E) při kolmém dopadu
r = (n1 - n2) / (n1 + n2);
t = 2*n1 / (n1 + n2);

% Energetické poměry (pro nemagnetická média, kolmý dopad)
R = abs(r)^2;
T = (n2/n1) * abs(t)^2;   % mělo by vycházet R+T = 1 (numericky s drobnou chybou)

% Amplitudy složek
A_inc = A0;          % incidentní vlna (→ +z) v n1
A_ref = r * A0;      % odražená vlna (→ -z) v n1
A_trn = t * A0;      % přenesená vlna (→ +z) v n2

%% --- Prostorová osa: několik λ na obou stranách rozhraní --------------------
L1 = 3*lambda1;              % rozsah vlevo (n1)
L2 = 3*lambda2;              % rozsah vpravo (n2)
NzL = 600; NzR = 600;

zL = linspace(-L1, 0, NzL);  % z < 0 v n1
zR = linspace(0,  L2, NzR);  % z >= 0 v n2
z  = [zL, zR(2:end)];        % spojení bez duplikace z=0

%% --- Rozvržení okna: vlevo graf, vpravo info panel -------------------------
fig = figure('Color','w','Name','EM wave: n1 -> n2 (normal incidence)');
tlo = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% Levý panel: pole E(z,t)
ax = nexttile(tlo,1); hold(ax,'on'); grid(ax,'on'); grid(ax,'minor');
ax.XLim = [min(z) max(z)];
% Y-limit zvol podle očekávané max. amplitudy vlevo (incident + odraz)
yMax = 1.2*max([abs(A_inc)+abs(A_ref), abs(A_trn)]);
ax.YLim = [-yMax yMax];
xlabel(ax,'z (m)'); ylabel(ax,'E (arb. u.)');
title(ax,'Šíření přes rozhraní (kolmý dopad)');

% Vykreslení rozhraní z=0
xline(ax,0,'--k','LineWidth',1);
text(ax, -0.02*range(ax.XLim), 0.88*yMax, sprintf('n_1 = %.3g',n1), ...
    'HorizontalAlignment','right');
text(ax,  +0.02*range(ax.XLim), 0.88*yMax, sprintf('n_2 = %.3g',n2), ...
    'HorizontalAlignment','left');

% Inicializace křivky (budeme jen přepisovat YData)
y0L = (A_inc)*sin(omega*0 - k1*zL + phi) + (showReflection*A_ref)*sin(omega*0 + k1*zL + phi);
y0R = (A_trn)*sin(omega*0 - k2*zR + phi);
y0  = [y0L, y0R(2:end)];
hWave = plot(ax, z, y0, 'b', 'LineWidth', 2);

% Pravý panel: info box (text mimo graf)
axInfo = nexttile(tlo,2); axis(axInfo,'off');
infoFmt = [ ...
    "GLOBÁLNÍ PARAMETRY", ...
    "-------------------", ...
    "c         = %.6g", ...
    "f         = %.6g  (λ_0 = %.6g)", ...
    "", ...
    "PROSTŘEDÍ 1 (z<0)", ...
    "-------------------", ...
    "n_1       = %.6g", ...
    "v_1       = %.6g", ...
    "λ_1       = %.6g", ...
    "A_{inc}   = %.6g", ...
    "A_{ref}   = %.6g   (r = %.6g)", ...
    "", ...
    "PROSTŘEDÍ 2 (z≥0)", ...
    "-------------------", ...
    "n_2       = %.6g", ...
    "v_2       = %.6g", ...
    "λ_2       = %.6g", ...
    "A_{trn}   = %.6g   (t = %.6g)", ...
    "", ...
    "ENERGIE", ...
    "-------------------", ...
    "R = |r|^2 = %.6g", ...
    "T = (n_2/n_1)|t|^2 = %.6g", ...
    "R + T     = %.6g", ...
    "", ...
    "STAV", ...
    "-------------------", ...
    "t         = %.2f s" ...
];
hInfo = text(axInfo, 0.02, 0.98, "", 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Consolas','FontSize',11, ...
    'BackgroundColor','w','Margin',6,'EdgeColor',[0.7 0.7 0.7]);

set(hInfo,'String', sprintf(strjoin(infoFmt,newline), ...
    c, f, lambda0, ...
    n1, v1, lambda1, A_inc, A_ref, r, ...
    n2, v2, lambda2, A_trn, t, ...
    R, T, R+T, ...
    0));

%% --- Animace ----------------------------------------------------------------
for i = 1:n_steps
    t = i*dt;

    % Vlevo (n1): incidentní → +z, odražená → -z (podle volby showReflection)
    yL = (A_inc)*sin(omega*t - k1*zL + phi) ...
       + (showReflection*A_ref)*sin(omega*t + k1*zL + phi);

    % Vpravo (n2): přenesená → +z
    yR = (A_trn)*sin(omega*t - k2*zR + phi);

    % Spoj do jedné křivky
    y = [yL, yR(2:end)];
    set(hWave, 'YData', y);

    % Přepiš info panel (jen čas)
    set(hInfo,'String', sprintf(strjoin(infoFmt,newline), ...
        c, f, lambda0, ...
        n1, v1, lambda1, A_inc, A_ref, r, ...
        n2, v2, lambda2, A_trn, t, ...
        R, T, R+T, ...
        t));

    title(ax, sprintf('Šíření přes rozhraní — t = %.2f s', t));
    drawnow limitrate;
    pause(Tpause/1000);
end
