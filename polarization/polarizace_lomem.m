% Polarizace opakovaným lomem: P(m,n) = m / ( m + (2n/(n^2 - 1))^2 )
clear; close all; clc;

P = @(m,n) m ./ ( m + (2*n./(n.^2 - 1)).^2 );

m = 1:50;                 % pocet pruchodu/odrazu
n_list = [1.33, 1.5, 2.4];% voda, sklo, diamant

figure('Color','w','Position',[200 200 800 500]); hold on;
for n = n_list
    plot(m, P(m,n), 'LineWidth', 2, 'DisplayName', sprintf('n = %.2f', n));
end
grid on; box on;
xlabel('Počet odrazů m','Interpreter','latex');
ylabel('Stupeň polarizace $P$','Interpreter','latex');
title('Stupeň polarizace $P$ v závislosti na počtu odrazů $m$','Interpreter','latex');
legend('Location','southeast','Interpreter','latex');
ylim([0 1]);

% Ulozit pro vlozeni do LaTeXu
exportgraphics(gcf, 'obrazky/polarizace_P_vs_m.png', 'Resolution', 300);
