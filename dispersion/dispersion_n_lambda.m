% dispersion_n_lambda.m
% Disperze indexu lomu n(λ) z Lorentzova oscilátoru (schematicky)

clear; clc; close all;

c = 3e8;                            % rychlost světla [m/s]
lam = linspace(300e-9, 900e-9, 1200);
omega = 2*pi*c./lam;

lam0   = 380e-9;                    % rezonanční vlnová délka [m]
omega0 = 2*pi*c/lam0;
gamma  = omega0/25;                 % tlumení
f      = 0.25;                      % síla (osc. síla) – bezrozměrné

% komplexní susceptibilita a permitivita
chi  = f*(omega0.^2) ./ (omega0.^2 - omega.^2 - 1i*gamma.*omega);
epsr = 1 + chi;

n_complex = sqrt(epsr);
n = real(n_complex);

fObj = saGetFigure('USER', [10 6]);
plot(lam*1e9, n, 'LineWidth', 2); grid on;
xlabel('Vlnová délka \lambda [nm]');
ylabel('Reálná část indexu lomu n(\lambda)');
title('Disperze indexu lomu (schematicky)');

% zvýraznění rezonanční oblasti
hold on;
w_nm = 100;                               % šířka zvýrazněné oblasti [nm]
xLo = (lam0 - w_nm*1e-9/2)*1e9;
xHi = (lam0 + w_nm*1e-9/2)*1e9;
yl = ylim;
patch([xLo xHi xHi xLo],[yl(1) yl(1) yl(2) yl(2)], [0.7 0.85 1], ...
      'FaceAlpha',0.15,'EdgeColor','none');

xline(lam0*1e9,'--');
saSaveFig(fObj, '../../Obrazky', 'disperze-index-lomu', 'png');
