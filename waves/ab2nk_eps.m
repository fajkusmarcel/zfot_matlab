function props = ab2nk_eps(alpha, beta, f, mu_r)
% Převod (alpha, beta, f, mu_r) -> n, kappa, eps_r', eps_r'', lambda
% Konvence: E(z) = E0 * exp(-alpha z) * exp(j(omega t - beta z))
if nargin < 4, mu_r = 1; end
c     = 299792458;
omega = 2*pi*f;

n     = (beta * c) / omega;       % lomová složka komplexního indexu
kappa = (alpha * c) / omega;      % extinkční složka (ztráty)

% Pozn.: pro mu_r ~= 1 platí (n - j kappa)^2 = eps_r * mu_r
eps_r_complex = (n - 1j*kappa).^2 / mu_r;
eps_r_p  = real(eps_r_complex);
eps_r_pp = -imag(eps_r_complex);  % ε = ε' - j ε''

% Vlnová délka v prostředí
lambda = 2*pi / beta * 1e9;

props = struct('n',n, ...
               'kappa',kappa, ...
               'eps_r_p',eps_r_p, ...
               'eps_r_pp',eps_r_pp, ...
               'lambda',lambda);
end
