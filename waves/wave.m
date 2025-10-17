%% Vlna v konkretnim bode (z = 0)
A = 1;          % amplituda
f = 1;          % frekvence [Hz]
omega = 2*pi*f; % uhlova frekvence [rad/s]
T = 1/f;        % perioda [s]
phi = 0;        % počáteční fáze [rad]

% Casova osa: 0..2T, jemne vzorkovani
t = linspace(0, 2, 1000);

% Casovy prubeh v bodu z=0
y = A * sin(omega * t + phi);

figure;
plot(t, y, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('y(t)');
title('Harmonicky kmit v case');


%% Vlna v prostoru (t = 0)
A = 1;            % amplituda          [-]
lambda = 1;     % vlnova delka       [m]
k = 2*pi/lambda;  % vlnove cislo       [rad/m]
phi = 0;         % pocatecni faze     [rad]

% Prostorova osa: 0..2*lambda, jemne vzorkovani
z = linspace(0, 2, 1000);

% Prostorovy prubeh v case t = 0
y = A * sin(k * z + phi);

figure;
plot(z, y, 'LineWidth', 1.5); grid on;
xlabel('z [m]'); ylabel('y(z)');
title('Harmonicka vlna v prostoru (t = 0)');

%% Postupna harmonicka vlna y(z,t) = A*sin(k*z - omega*t + phi)
A = 1;              % amplituda
lambda = 1.0;       % vlnova delka [m]
k = 2*pi/lambda;    % vlnove cislo [rad/m]
f = 1.0;            % frekvence [Hz]
omega = 2*pi*f;     % uhlova frekvence [rad/s]
phi = 0;           % pocatecni faze [rad]

% Prostor a cas pro "snapshoty"
z = linspace(0, 3, 1000);
t_list = [0, 0.10, 0.20, 0.30];    % casy v sekundach (zmente a sledujte posun)

figure; hold on; grid on
for t = t_list
    y = A * sin(k*z - omega*t + phi);
    plot(z, y, 'LineWidth', 1.3);
end

y = A * sin(k*z - omega*t_list(1) + phi);  plot(z, y, 'LineWidth', 1.3);
y = A * sin(k*z - omega*t_list(2) + phi);  plot(z, y, 'LineWidth', 1.3);
y = A * sin(k*z - omega*t_list(3) + phi);  plot(z, y, 'LineWidth', 1.3);
y = A * sin(k*z - omega*t_list(4) + phi);  plot(z, y, 'LineWidth', 1.3);

xlabel('z [m]'); ylabel('y(z,t)');
title('Postupna vlna: profily v ruznych casech (smer +z)');
legend(arrayfun(@(tt) sprintf('t = %.2f s', tt), t_list, 'UniformOutput', false), ...
       'Location', 'southoutside');
hold off


%% Postupna harmonicka vlna - simulace y(z,t) = A*sin(k*z - omega*t + phi0)
A = 1;              % amplituda
lambda = 1.0;       % vlnova delka [m]
k = 2*pi/lambda;    % vlnove cislo [rad/m]
f = 1.0;            % frekvence [Hz]
omega = 2*pi*f;     % uhlova frekvence [rad/s]
phi0 = 0;           % pocatecni faze [rad]

% Prostor a cas pro "snapshoty"
z = linspace(0, 3, 100);
t_list = linspace(0,1,100) % casy v sekundach (zmente a sledujte posun)

figure; hold off; grid on
for t = t_list
    y = A * sin(k*z - omega*t + phi0);
    plot(z, y, 'LineWidth', 1.3);
    pause(0.1)
end

xlabel('z [m]'); ylabel('y(z,t)');
title('Postupna vlna: profily v ruznych casech (smer +z)');
legend(arrayfun(@(tt) sprintf('t = %.2f s', tt), t_list, 'UniformOutput', false), ...
       'Location', 'southoutside');
hold off


