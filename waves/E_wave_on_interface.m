%% Harmonická vlna ve dvou prostředích + výpočet n, kappa, eps'r, eps''r
clear; clc; close all;

% ===== Parametry =====
E0   = 1;               % amplituda [V/m]
f    = 2e14;            % frekvence [Hz]  (např. ~200 THz ~ 1.5 um; uprav dle potřeby)
t    = 0;               % časový okamžik pro snapshot
mu_r = 1;               % relativní permeabilita (většina dielektrik ~1)

% Prostředí 1 (vlevo, z ∈ [-N*lambda1, 0]) – bezeztrátové
beta1  = 1;             % [rad/m]
alpha1 = 0;             % [1/m]
nPerL  = 3;             % kolik period vlevo zobrazit

% Prostředí 2 (vpravo, z ∈ [0, N*lambda2]) – ztrátové
beta2  = 3;             % [rad/m]
alpha2 = 0.25;          % [1/m] mírnější útlum pro hezčí průběh
nPerR  = 6;             % kolik period vpravo zobrazit

% ===== Odvozené veličiny =====
c     = 299792458;
omega = 2*pi*f;

lambda1 = 2*pi/beta1;       % zadané beta -> "geometrická" λ1
lambda2 = 2*pi/beta2;       % zadané beta -> "geometrická" λ2

% Vypočti n, kappa, eps_r', eps_r'' z (alpha, beta, f)
p1 = ab2nk_eps(alpha1, beta1, f, mu_r);
p2 = ab2nk_eps(alpha2, beta2, f, mu_r);

% ===== Z-osa a pole =====
z1 = linspace(-nPerL*lambda1, 0, 4000);
E1 = E0 .* exp(-alpha1.*(z1 - z1(end))) .* cos(omega*t - beta1.*z1);  % alpha1=0

z2 = linspace(0, nPerR*lambda2, 8000);
E2 = E0 .* exp(-alpha2.*z2) .* cos(omega*t - beta2.*z2);

% ===== Vykreslení =====
figure('Color','w'); hold on;

% podbarvení pravé oblasti
yspan = 1.25*E0;
patch([0 max(z2) max(z2) 0], [-yspan -yspan yspan yspan], ...
      [0.93 0.93 0.93], 'EdgeColor','none', 'HandleVisibility','off');

plot(z1, E1, 'b', 'LineWidth', 1.2);
plot(z2, E2, 'b', 'LineWidth', 1.2);

% dělicí čára
plot([0 0], [-yspan yspan], 'k', 'LineWidth', 1);

xlim([min(z1) max(z2)]);
ylim([-yspan yspan]);
grid on; box on;
xlabel('z (m)');
ylabel('Field Amplitude (V/m)');
title('Harmonická vlna: bezeztrátové (vlevo) a ztrátové (vpravo) prostředí');

% Anotace parametrů (beta, alpha, n, kappa) – s LaTeX interpreterem
txtL = sprintf('$\\beta_1 = %.3g\\,\\mathrm{rad/m},\\; \\alpha_1 = %.3g\\,\\mathrm{m^{-1}}$ \n $n_1 = %.3g,\\; \\kappa_1 = %.3g$', ...
               beta1, alpha1, p1.n, p1.kappa);

txtR = sprintf('$\\beta_2 = %.3g\\,\\mathrm{rad/m},\\; \\alpha_2 = %.3g\\,\\mathrm{m^{-1}}$ \n $n_2 = %.3g,\\; \\kappa_2 = %.3g$', ...
               beta2, alpha2, p2.n, p2.kappa);

% Vykreslení textu do grafu
text(min(z1)+0.02*range([min(z1) max(z2)]),  0.88*yspan, txtL, ...
    'FontSize', 10, 'Interpreter','latex');

text(0.02*max(z2), 0.88*yspan, txtR, ...
    'FontSize', 10, 'Interpreter','latex');


% ===== Výpis do konzole =====
disp('--- Prostředí 1 ---');
disp(p1);
disp('--- Prostředí 2 ---');
disp(p2);

%% ===== Pomocná funkce =====

