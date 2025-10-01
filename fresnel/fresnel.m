% Fresnelovy koeficienty - porovnani obecneho a trigonometrickeho tvaru

%% ========================================================================
% Fresnelovy relace pro n_1 < n_2
close all;, clc; clear all;

% Indexy lomu
n1 = 1.0;     % např. vzduch
n2 = 1.5;     % např. sklo


% Rozsah úhlu dopadu (0 až těsně pod 90°)
theta_i = linspace(0, pi/2-0.001, 20);

% Úhel lomu z Snellova zákona
theta_t = asin(n1./n2 .* sin(theta_i));

% Brewsterův úhel
theta_B = atan(n2/n1);   % [rad]



% --- Obecné tvary Fresnelových vztahů ---
rs = (n1*cos(theta_i) - n2*cos(theta_t)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
rp = (n2*cos(theta_i) - n1*cos(theta_t)) ./ (n2*cos(theta_i) + n1*cos(theta_t));

ts = (2*n1*cos(theta_i)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
tp = (2*n1*cos(theta_i)) ./ (n2*cos(theta_i) + n1*cos(theta_t));

fObj = saGetFigure('USER', [18 8]);

subplot(1,2,1);

plot(theta_i*180/pi, rs, 'r-', 'LineWidth', 1.5); hold on;
plot(theta_i*180/pi, rp, 'b-', 'LineWidth', 1.5); hold on;

% vodorovná čára v nule
yline(0, 'k--', 'LineWidth', 1);


% svislá čára v Brewsterově úhlu
xline(theta_B*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_B','LabelOrientation','horizontal');


xlabel('\theta_i [°]');
ylabel('Reflexní koeficient');
legend('r_s','r_p','Location','southwest');
grid minor;
xlim([0 90])
title('Reflexní koeficenty');

subplot(1,2,2);
plot(theta_i*180/pi, ts, 'r-', 'LineWidth', 1.5); hold on;
plot(theta_i*180/pi, tp, 'b-', 'LineWidth', 1.5); hold on;
xlabel('\theta_i [°]');
ylabel('Transmisní koeficient');
legend('t_s','t_p','Location','southwest');
grid minor;
xlim([0 90])
title('Transmisní koeficienty');

saSaveFig(fObj, '../../Obrazky/', 'fresnel-koeficienty-a', 'png');

%% ========================================================================
% Fresnelovy relace pro n_1 > n_2
close all;, clc; clear all;

% Indexy lomu
n1 = 1.5;     % např. vzduch
n2 = 1.0;     % např. sklo


% Rozsah úhlu dopadu (0 až těsně pod 90°)
theta_i = linspace(0, pi/2-0.001, 2000);

% Úhel lomu z Snellova zákona
theta_t = asin(n1./n2 .* sin(theta_i));

% Brewsterův úhel
theta_B = atan(n2/n1);   % [rad]
% Kriticky úhel
theta_C = asin(n2/n1);   % [rad]


% --- Obecné tvary Fresnelových vztahů ---
rs = (n1*cos(theta_i) - n2*cos(theta_t)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
rp = (n2*cos(theta_i) - n1*cos(theta_t)) ./ (n2*cos(theta_i) + n1*cos(theta_t));

ts = (2*n1*cos(theta_i)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
tp = (2*n1*cos(theta_i)) ./ (n2*cos(theta_i) + n1*cos(theta_t));

fObj = saGetFigure('USER', [18 8]);

subplot(1,2,1);

mask = theta_i <= theta_C;
rs = rs(mask);
rp = rp(mask);
ts = rs(mask);
tp = rp(mask);
theta_i = theta_i(mask);

plot(theta_i*180/pi, rs, 'r-', 'LineWidth', 1.5); hold on;
plot(theta_i*180/pi, rp, 'b-', 'LineWidth', 1.5); hold on;

% vodorovná čára v nule
yline(0, 'k--', 'LineWidth', 1);


% svislá čára v Brewsterově úhlu
xline(theta_B*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_B','LabelOrientation','horizontal');
% svislá čára v kritickem úhlu
xline(theta_C*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_C','LabelOrientation','horizontal');

xlabel('\theta_i [°]');
ylabel('Reflexní koeficient');
legend('r_s','r_p','Location','southwest');
grid minor;
xlim([0 90])
title('Reflexní koeficenty');

subplot(1,2,2);
plot(theta_i*180/pi, ts, 'r-', 'LineWidth', 1.5); hold on;
plot(theta_i*180/pi, tp, 'b-', 'LineWidth', 1.5); hold on;
xlabel('\theta_i [°]');
ylabel('Transmisní koeficient');
legend('t_s','t_p','Location','southwest');
grid minor;
xlim([0 90])
title('Transmisní koeficienty');

saSaveFig(fObj, '../../Obrazky/', 'fresnel-koeficienty-b', 'png');




%% ========================================================================
% Výkonové Fresnelovy koeficienty pro n_1 < n_2
close all; clc; clear;

% Indexy lomu
n1 = 1.0;      % vzduch
n2 = 1.5;      % sklo

% Úhly
theta_i = linspace(0, pi/2-1e-3, 500);              % dopad
theta_t = asin(n1./n2 .* sin(theta_i));             % lom (Snell)
theta_B = atan(n2/n1);                               % Brewster [rad]

% Amplitudové koeficienty (stejné jako ve tvém skriptu)
rs = (n1*cos(theta_i) - n2*cos(theta_t)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
rp = (n2*cos(theta_i) - n1*cos(theta_t)) ./ (n2*cos(theta_i) + n1*cos(theta_t));
ts = (2*n1*cos(theta_i)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
tp = (2*n1*cos(theta_i)) ./ (n2*cos(theta_i) + n1*cos(theta_t));

% Výkonové (intenzitové) koeficienty
Rs = abs(rs).^2;
Rp = abs(rp).^2;
Ts = (n2.*cos(theta_t)./(n1.*cos(theta_i))) .* abs(ts).^2;
Tp = (n2.*cos(theta_t)./(n1.*cos(theta_i))) .* abs(tp).^2;

% (volitelné) kontrola zákona zachování energie
% max(abs([Rs+Ts, Rp+Tp] - 1), [], 'all')

% --- Kreslení ---
fObj = saGetFigure('USER', [18 8]);

% Odraz
subplot(1,2,1);
plot(theta_i*180/pi, Rs, 'r-', 'LineWidth', 1.6); hold on;
plot(theta_i*180/pi, Rp, 'b-', 'LineWidth', 1.6);
yline(1, 'k:', 'LineWidth', 1);           % referenční 100 %
xline(theta_B*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_B','LabelOrientation','horizontal');
xlabel('\theta_i [°]'); ylabel('Výkonový koeficient odrazu');
legend('R_s','R_p','Location','northwest');
title('Výkonové odrazové koeficienty');
grid minor; xlim([0 90]); ylim([0 1]);

% Přenos
subplot(1,2,2);
plot(theta_i*180/pi, Ts, 'r-', 'LineWidth', 1.6); hold on;
plot(theta_i*180/pi, Tp, 'b-', 'LineWidth', 1.6);
xline(theta_B*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_B','LabelOrientation','horizontal');
xlabel('\theta_i [°]'); ylabel('Výkonový koeficient přenosu');
legend('T_s','T_p','Location','northeast');
title('Výkonové transmisní koeficienty');
grid minor; xlim([0 90]); ylim([0 1]);

saSaveFig(fObj, '../../Obrazky/', 'fresnel-vykon-a', 'png');


%% ========================================================================
% Výkonové Fresnelovy koeficienty pro n_1 > n_2
close all; clc; clear;

% Indexy lomu
n1 = 1.5;      % vzduch
n2 = 1.0;      % sklo

% Úhly
theta_i = linspace(0, pi/2-1e-3, 500);              % dopad
theta_t = asin(n1./n2 .* sin(theta_i));             % lom (Snell)
theta_B = atan(n2/n1);                               % Brewster [rad]
% Kriticky úhel
theta_C = asin(n2/n1);   % [rad]

% Amplitudové koeficienty (stejné jako ve tvém skriptu)
rs = (n1*cos(theta_i) - n2*cos(theta_t)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
rp = (n2*cos(theta_i) - n1*cos(theta_t)) ./ (n2*cos(theta_i) + n1*cos(theta_t));
ts = (2*n1*cos(theta_i)) ./ (n1*cos(theta_i) + n2*cos(theta_t));
tp = (2*n1*cos(theta_i)) ./ (n2*cos(theta_i) + n1*cos(theta_t));

% Výkonové (intenzitové) koeficienty
Rs = abs(rs).^2;
Rp = abs(rp).^2;
Ts = (n2.*cos(theta_t)./(n1.*cos(theta_i))) .* abs(ts).^2;
Tp = (n2.*cos(theta_t)./(n1.*cos(theta_i))) .* abs(tp).^2;

% Maskování

mask = theta_i <= theta_C + 0.01;
Rs = Rs(mask);
Rp = Rp(mask);
Ts = Ts(mask);
Tp = Tp(mask);
theta_i = theta_i(mask);

% (volitelné) kontrola zákona zachování energie
% max(abs([Rs+Ts, Rp+Tp] - 1), [], 'all')

% --- Kreslení ---
fObj = saGetFigure('USER', [18 8]);

% Odraz
subplot(1,2,1);
plot(theta_i*180/pi, Rs, 'r-', 'LineWidth', 1.6); hold on;
plot(theta_i*180/pi, Rp, 'b-', 'LineWidth', 1.6);
yline(1, 'k:', 'LineWidth', 1);           % referenční 100 %

% svislá čára v kritickem úhlu
xline(theta_C*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_C','LabelOrientation','horizontal');

xline(theta_B*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_B','LabelOrientation','horizontal');
xlabel('\theta_i [°]'); ylabel('Výkonový koeficient odrazu');
legend('R_s','R_p','Location','northwest');
title('Výkonové odrazové koeficienty');
grid minor; xlim([0 90]); ylim([0 1]);

% Přenos
subplot(1,2,2);
plot(theta_i*180/pi, Ts, 'r-', 'LineWidth', 1.6); hold on;
plot(theta_i*180/pi, Tp, 'b-', 'LineWidth', 1.6);
xline(theta_B*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_B','LabelOrientation','horizontal');
% svislá čára v kritickem úhlu
xline(theta_C*180/pi, 'k--', 'LineWidth', 1.2, ...
      'Label','\theta_C','LabelOrientation','horizontal');
xlabel('\theta_i [°]'); ylabel('Výkonový koeficient přenosu');
legend('T_s','T_p','Location','northeast');
title('Výkonové transmisní koeficienty');
grid minor; xlim([0 90]); ylim([0 1]);

saSaveFig(fObj, '../../Obrazky/', 'fresnel-vykon-b', 'png');

