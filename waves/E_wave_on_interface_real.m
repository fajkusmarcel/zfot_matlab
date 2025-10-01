%% E_wave_on_interface.m
% Rovinná harmonická vlna na rozhraní dvou prostředí
clear; clc; close all;

% ===== Zadání =====
E0   = 1;               % amplituda [V/m]
f    = 200e12;          % frekvence [Hz] (200 THz ~ 1.5 µm)
t    = 0;               % časový snapshot
mu_r = 1;               % relativní permeabilita
% =========================================================
% Komplexní index lomu:  ñ = n - jκ
%
% - Reálná část (n):    určuje fázovou rychlost a tím vlnovou délku
%                       λ = λ0 / n, kde λ0 = c / f.
%                       Čím větší n, tím kratší vlnová délka v prostředí.
%
% - Imaginární část (κ): tzv. extinkční koeficient.
%                       Určuje absorpci (tlumení) vlny:
%                       α = (ω/c)·κ  [m⁻¹].
%                       Čím větší κ, tím rychlejší exponenciální pokles amplitudy.
%
% =========================================================

% Prostředí 1
n1     = 1;             % Reálná část indexu lomu (odpovídá vakuu/luftu).
                        % Vlna se šíří rychlostí c, vlnová délka λ = λ0.
kappa1 = 1e-6;          % Imaginární část indexu lomu (prakticky nulová).
                        % Tlumení je zanedbatelné → α ≈ 0, vlna se netlumí.

% Prostředí 2
n2     = 1.5;           % Reálná část indexu lomu.
                        % Vlna se šíří poloviční rychlostí oproti vakuu,
                        % vlnová délka se zkrátí na λ = λ0 / 2.
kappa2 = 0.05;          % Imaginární část indexu lomu.
                        % Útlum je výrazný, vlna rychle zaniká při průchodu prostředím.
                        % Příklad: absorbující dielektrikum, polovodič, barvené sklo.


% ===== Výpočet parametrů =====
c     = 299792458;
omega = 2*pi*f;

% Prostředí 1
beta1  = omega*n1/c;
alpha1 = omega*kappa1/c;
lambda1 = 2*pi/beta1;

% Prostředí 2
beta2  = omega*n2/c;
alpha2 = omega*kappa2/c;
lambda2 = 2*pi/beta2;

% ===== Pole =====
nPerL = 3; 
nPerR = 6; 

z1 = linspace(-nPerL*lambda1, 0, 2000);
E1 = E0 .* exp(-alpha1.*(z1 - z1(end))) .* cos(omega*t - beta1.*z1);

z2 = linspace(0, nPerR*lambda2, 4000);
E2 = E0 .* exp(-alpha2.*z2) .* cos(omega*t - beta2.*z2);

% ===== Vykreslení =====
figure('Color','w'); hold on;
yspan = 1.5*E0;

patch([0 max(z2) max(z2) 0], [-yspan -yspan yspan yspan], ...
      [0.93 0.93 0.93], 'EdgeColor','none', 'HandleVisibility','off');

plot(z1, E1, 'b', 'LineWidth', 1.2);
plot(z2, E2, 'b', 'LineWidth', 1.2);

plot([0 0], [-yspan yspan], 'k', 'LineWidth', 1);

xlim([min(z1) max(z2)]);
ylim([-yspan yspan]);
grid on; box on;
xlabel('z (m)');
ylabel('Field Amplitude (V/m)');
title('Harmonická vlna: prostředí 1 a 2');

% ===== Anotace =====
txtL = sprintf('$n_1 = %.3f$ \n $\\kappa_1 = %.1e$ \n $\\beta_1 = %.2e\\,\\mathrm{rad/m}$ \n $\\alpha_1 = %.2e\\,\\mathrm{m^{-1}}$', ...
               n1, kappa1, beta1, alpha1);

txtR = sprintf('$n_2 = %.3f$ \n $\\kappa_2 = %.3f$ \n $\\beta_2 = %.2e\\,\\mathrm{rad/m}$ \n $\\alpha_2 = %.2e\\,\\mathrm{m^{-1}}$', ...
               n2, kappa2, beta2, alpha2);

text(min(z1)+0.02*range([min(z1) max(z2)]),  0.8*yspan, txtL, ...
    'FontSize', 10, 'Interpreter','latex');

text(0.02*max(z2), 0.8*yspan, txtR, ...
    'FontSize', 10, 'Interpreter','latex');


% Vypočti n, kappa, eps_r', eps_r'' z (alpha, beta, f)
p1 = ab2nk_eps(alpha1, beta1, f, mu_r);
p2 = ab2nk_eps(alpha2, beta2, f, mu_r);
% ===== Výpis do konzole =====
disp('--- Prostředí 1 ---');
disp(p1);
disp('--- Prostředí 2 ---');
disp(p2);
