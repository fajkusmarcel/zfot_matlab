%% Animace šíření elektrické složky EM vlny v 1D prostoru – info panel mimo graf

clear; clc; close all;

%% Parametry
c = 1;                           % rychlost světla ve vakuu (normalizovaně)
z = linspace(-1, 1, 800);        % prostorová osa [m]
n = 1.5;                         % index lomu
v = c / n;                       % rychlost šíření
f   = 1.00;                      % frekvence [Hz]
phi = 0.00;                      % počáteční fáze [rad]
A   = 1.00;                      % amplituda

lambda = v / f;                  % vlnová délka
omega  = 2*pi*f;                 % úhlová frekvence
k      = 2*pi / lambda;          % vlnové číslo

n_steps = 200;                   % kroky animace
dt      = 0.01;                  % časový krok [s]
Tpause  = 50;                    % pauza [ms]
z_point = 0;                     % sledovaný bod

%% Rozvržení: vlevo graf, vpravo info panel
fig = figure('Color','w','Name','EM wave animation');
tlo = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% --- Levý panel: graf ---
ax = nexttile(tlo,1); hold(ax,'on'); grid(ax,'on'); grid(ax,'minor');
ax.XLim = [min(z) max(z)];
ax.YLim = [-1.1*A 1.1*A];
xlabel(ax,'z (m)'); ylabel(ax,'E (arb. u.)');
title(ax,'Animace postupné sinusové vlny');

y0 = A * sin(omega*0 - k*z + phi);
hWave  = plot(ax, z, y0, 'b', 'LineWidth', 2);
yP0    = A * sin(omega*0 - k*z_point + phi);
hPoint = plot(ax, z_point, yP0, 'r.', 'MarkerSize', 36);

% --- Pravý panel: čistě text (osy vypnuté) ---
axInfo = nexttile(tlo,2); axis(axInfo,'off');
% Jeden textový objekt, který budeme přepisovat (monospace font pro hezké zarovnání)
hInfo = text(axInfo, 0.02, 0.98, "", 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Consolas','FontSize',11, ...
    'BackgroundColor','w','Margin',6,'EdgeColor',[0.7 0.7 0.7]);

% Úvodní obsah panelu
set(hInfo,'String', sprintf([ ...
    'PARAMETRY\n' ...
    '-----------\n' ...
    'n        = %.3g\n' ...
    'v        = %.3g m/s\n' ...
    'f        = %.3g Hz\n' ...
    'A        = %.3g\n' ...
    'λ        = %.3g m\n' ...
    '\nSTAV\n' ...
    '-----------\n' ...
    't        = %.2f s\n' ...
    'φ(0,t)   = %.2f rad (%.1f°)\n' ...
], n, v, f, A, lambda, 0, 0, 0));

%% Smyčka animace
for i = 1:n_steps
    t = i*dt;

    % Vlna a bod v z=0
    y = A * sin(omega*t - k*z + phi);
    y_point = A * sin(omega*t - k*z_point + phi);

    % Update křivek
    set(hWave,  'YData', y);
    set(hPoint, 'YData', y_point);

    % Aktuální fáze v z=0 (zabalená do [0, 2π) jen pro výpis)
    phi0 = mod(omega*t - k*z_point + phi, 2*pi);

    % Přepiš info panel (vpravo)
    set(hInfo,'String', sprintf([ ...
        'PARAMETRY\n' ...
        '-----------\n' ...
        'n        = %.3g\n' ...
        'v        = %.3g m/s\n' ...
        'f        = %.3g Hz\n' ...
        'A        = %.3g\n' ...
        'λ        = %.3g m\n' ...
        '\nSTAV\n' ...
        '-----------\n' ...
        't        = %.2f s\n' ...
        'φ(0,t)   = %.2f rad (%.1f°)\n' ...
    ], n, v, f, A, lambda, t, phi0, phi0*180/pi));

    title(ax, sprintf('Animace postupné sinusové vlny  —  t = %.2f s', t));

    drawnow limitrate;
    pause(Tpause/1000);
end
