%% 
close all;
alpha = 0:0.01:0.8;


y_pi = alpha;
y_sin = sin(alpha);
y_tan = tan(alpha);



% Nastavení grafiky
fObj = saGetFigure('USER', [10 6]);
fObj.xlabel = 'Úhel (rad)';
fObj.ylabel = 'Hodnota (-)';
fObj.legend = {'$\alpha$', '$\sin \alpha$', '$\tan \alpha$', 'Odchylka'};
fObj.Interpreter = 'latex';
fObj.FontSize = 10;
fObj.legendFontSize = 8;
fObj.legendLocation = 'northwest';

plot(alpha, y_pi); hold on;
plot(alpha, y_sin);
plot(alpha, y_tan);
plot(alpha, y_pi - y_sin);

plot([0.3 0.3], ylim, '--k'); % od x=0.3 po celý rozsah osy y

xline(0.3, '--k', '$0.3\,\mathrm{rad} \approx 17^\circ$', ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'top', ...
    'Interpreter', 'latex');



saSetFigure(fObj);
saSaveFig(fObj, '../../Obrazky', 'paraxial', 'png');




