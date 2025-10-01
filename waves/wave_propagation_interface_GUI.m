%% EM vlna přes rozhraní n1 -> n2 (kolmý dopad) s posuvníky
% Vlevo: uiaxes s animací (incident + volitelně odražená) v n1 a přenesená v n2
% Vpravo: infopanel a ovládací prvky (f, n1, n2, A0, Show reflection)

clear; clc; close all;

%% --- Výchozí parametry ---
S.c   = 1;          % rychlost světla ve vakuu (normalizováno)
S.f   = 1.00;       % frekvence [Hz] (společná ve všech prostředích)
S.A0  = 1.00;       % amplituda incidentní vlny
S.phi = 0.00;       % počáteční fáze [rad]
S.n1  = 1.00;       % index lomu vlevo (z<0)
S.n2  = 1.50;       % index lomu vpravo (z>=0)
S.showReflection = false;
S.showFresnell = false;
S.useAttn2 = false;   % zap/vyp útlum v n2
S.alpha2   = 0.25;    % koeficient útlumu amplitudy v n2 [Np/m]

% Animace
S.dt     = 0.01;    % časový krok [s]
S.Tpause = 40;      % pauza mezi snímky [ms]

%% --- UI rozložení --- (puvodni verze)
% uif = uifigure('Name','EM wave: n1 -> n2 (normal incidence)','Color','w','Position',[100 100 1100 800]);
% ax  = uiaxes(uif,'Position',[40 90 640 650]); hold(ax,'on'); grid(ax,'on'); ax.GridAlpha = 0.25;
% xlabel(ax,'z (m)'); ylabel(ax,'E (arb. u.)'); title(ax,'Šíření přes rozhraní (kolmý dopad)');
% 
% % Infopanel (needitovatelný)
% ta = uitextarea(uif,'Position',[710 320 320 420], 'Editable','off', ...
%     'FontName','Consolas','FontSize',11);
% 
% % Ovládání: štítky + posuvníky/checkbox
% % Frekvence
% uilabel(uif,'Text','f [Hz]','Position',[710 290 60 22]);
% slF = uislider(uif,'Position',[770 310 260 3],'Limits',[0.10 5.00],'Value',S.f);
% lblF = uilabel(uif,'Text',sprintf('%.2f',S.f),'Position',[1035 290 50 22]);
% 
% % n1
% uilabel(uif,'Text','n_1','Position',[710 250 60 22]);
% slN1 = uislider(uif,'Position',[770 260 260 3],'Limits',[1.00 2.50],'Value',S.n1);
% lblN1 = uilabel(uif,'Text',sprintf('%.3f',S.n1),'Position',[1035 250 50 22]);
% 
% % n2
% uilabel(uif,'Text','n_2','Position',[710 210 60 22]);
% slN2 = uislider(uif,'Position',[770 220 260 3],'Limits',[1.00 2.50],'Value',S.n2);
% lblN2 = uilabel(uif,'Text',sprintf('%.3f',S.n2),'Position',[1035 210 50 22]);
% 
% % Amplituda
% uilabel(uif,'Text','A_0','Position',[710 170 60 22]);
% slA = uislider(uif,'Position',[770 180 260 3],'Limits',[0.10 2.00],'Value',S.A0);
% lblA = uilabel(uif,'Text',sprintf('%.2f',S.A0),'Position',[1035 170 50 22]);
% 
% % Label + slider + numeric input pro alpha2
% uilabel(uif,'Text','\alpha_2 [Np/m]','Position',[710 160 80 22]);
% slAlpha2 = uislider(uif,'Position',[710 150 320 3], 'Limits',[0 2], 'Value',S.alpha2);
% S.slAlpha2.ValueChangingFcn = @(src,evt) onChange('ALPHA2', evt.Value);
% 
% efAlpha2 = uieditfield(uif,'numeric', 'Position',[1035 160 70 22], 'Limits',[0 Inf], 'Value',S.alpha2);
% S.efAlpha2.ValueChangedFcn = @(src,evt) onChange('ALPHA2', evt.Value);
% 
% % Checkbox odrazu
% cbF = uicheckbox(uif,'Text','Zahrnout Fresnellovy rovnice','Position',[710 100 200 24],'Value',S.showFresnell);
% cbR = uicheckbox(uif,'Text','Zobrazit odraženou vlnu','Position',[710 70 200 24],'Value',S.showReflection);
% cbAttn2 = uicheckbox(uif,'Text','Zohlednit útlum v n_2', 'Position',[710 40 220 24], 'Value',S.useAttn2);
% 
% Ulož GUI handlery a grafické objekty do S
% S.uif = uif; S.ax = ax; S.ta = ta;
% S.slF = slF; S.slN1 = slN1; S.slN2 = slN2; S.slA = slA;
% S.lblF = lblF; S.lblN1 = lblN1; S.lblN2 = lblN2; S.lblA = lblA; S.cbF = cbF; S.cbR = cbR;S.cbAttn2 = cbAttn2;S.slAlpha2 = slAlpha2;S.efAlpha2 = efAlpha2;



%% --- UI rozložení (responsivní přes uigridlayout) ---

uif = uifigure('Name','EM wave: n1 -> n2 (normal incidence)', ...
               'Color','w','Position',[100 100 1100 800]);

% Hlavní grid: 1 řádek, 2 sloupce (vlevo graf | vpravo panel)
g = uigridlayout(uif,[1 2]);
g.ColumnWidth   = {'3x','1.2x'};
g.RowHeight     = {'1x'};
g.Padding       = [16 16 16 16];
g.ColumnSpacing = 16;

% --- Levý sloupec: osy ---
ax = uiaxes(g); hold(ax,'on'); grid(ax,'on'); ax.GridAlpha = 0.25;
ax.Layout.Row = 1; ax.Layout.Column = 1;
xlabel(ax,'z (m)'); ylabel(ax,'E (arb. u.)');
title(ax,'Šíření přes rozhraní (kolmý dopad)');

% --- Pravý sloupec: panel s ovládáním a info ---
p = uipanel(g,'Title','Parametry a info');
p.Layout.Row = 1; p.Layout.Column = 2;

% Vnitřní grid v panelu (2 sloupce pro label|hodnota a řádky pro prvky)
pg = uigridlayout(p,[13 2]);
pg.ColumnWidth = {'1x','fit'};
pg.RowHeight = {'fit','fit', ...  % f (label+value; slider)
                'fit','fit', ...  % n1
                'fit','fit', ...  % n2
                'fit','fit', ...  % A0
                'fit','fit', ...  % alpha2 (label+edit; slider)
                'fit', ...        % checkboxes řádek 1 (Fresnel + Reflection)
                'fit', ...        % checkbox řádek 2 (Attn n2)
                '1x'};            % textarea zabere zbytek

% === f [Hz] ===
lblFName = uilabel(pg,'Text','f [Hz]','FontWeight','bold');
lblFName.Layout.Row = 1; lblFName.Layout.Column = 1;
lblF = uilabel(pg,'Text',sprintf('%.2f',S.f));
lblF.Layout.Row = 1; lblF.Layout.Column = 2;

slF = uislider(pg,'Limits',[0.10 5.00],'Value',S.f);
slF.Layout.Row = 2; slF.Layout.Column = [1 2];
slF.ValueChangingFcn = @(src,evt) onChange('f',evt.Value);

% === n1 ===
lblN1Name = uilabel(pg,'Text','n_1','FontWeight','bold');
lblN1Name.Layout.Row = 3; lblN1Name.Layout.Column = 1;
lblN1 = uilabel(pg,'Text',sprintf('%.3f',S.n1));
lblN1.Layout.Row = 3; lblN1.Layout.Column = 2;

slN1 = uislider(pg,'Limits',[1.00 2.50],'Value',S.n1);
slN1.Layout.Row = 4; slN1.Layout.Column = [1 2];
slN1.ValueChangingFcn = @(src,evt) onChange('n1',evt.Value);

% === n2 ===
lblN2Name = uilabel(pg,'Text','n_2','FontWeight','bold');
lblN2Name.Layout.Row = 5; lblN2Name.Layout.Column = 1;
lblN2 = uilabel(pg,'Text',sprintf('%.3f',S.n2));
lblN2.Layout.Row = 5; lblN2.Layout.Column = 2;

slN2 = uislider(pg,'Limits',[1.00 2.50],'Value',S.n2);
slN2.Layout.Row = 6; slN2.Layout.Column = [1 2];
slN2.ValueChangingFcn = @(src,evt) onChange('n2',evt.Value);

% === A0 ===
lblAName = uilabel(pg,'Text','A_0','FontWeight','bold');
lblAName.Layout.Row = 7; lblAName.Layout.Column = 1;
lblA = uilabel(pg,'Text',sprintf('%.2f',S.A0));
lblA.Layout.Row = 7; lblA.Layout.Column = 2;

slA = uislider(pg,'Limits',[0.10 2.00],'Value',S.A0);
slA.Layout.Row = 8; slA.Layout.Column = [1 2];
slA.ValueChangingFcn = @(src,evt) onChange('A0',evt.Value);

% === alpha2 ===
lblAlpha2Name = uilabel(pg,'Text','α₂ [Np/m]','FontWeight','bold'); % klidně 'alpha_2'
lblAlpha2Name.Layout.Row = 9; lblAlpha2Name.Layout.Column = 1;

efAlpha2 = uieditfield(pg,'numeric','Limits',[0 Inf],'Value',S.alpha2);
efAlpha2.Layout.Row = 9; efAlpha2.Layout.Column = 2;
efAlpha2.ValueChangedFcn = @(src,evt) onChange('ALPHA2',evt.Value);

slAlpha2 = uislider(pg,'Limits',[0 2],'Value',S.alpha2);
slAlpha2.Layout.Row = 10; slAlpha2.Layout.Column = [1 2];
slAlpha2.ValueChangingFcn = @(src,evt) onChange('ALPHA2',evt.Value);

% === Checkboxes (2 řádky) ===
cbF = uicheckbox(pg,'Text','Zahrnout Fresnelovy koeficienty','Value',S.showFresnell);
cbF.Layout.Row = 11; cbF.Layout.Column = 1;
cbF.ValueChangedFcn = @(src,evt) onChange('F',evt.Value);

cbR = uicheckbox(pg,'Text','Zobrazit odraženou vlnu','Value',S.showReflection);
cbR.Layout.Row = 11; cbR.Layout.Column = 2;
cbR.ValueChangedFcn = @(src,evt) onChange('R',evt.Value);

cbAttn2 = uicheckbox(pg,'Text','Zohlednit útlum v n_2','Value',S.useAttn2);
cbAttn2.Layout.Row = 12; cbAttn2.Layout.Column = [1 2];
cbAttn2.ValueChangedFcn = @(src,evt) onChange('ATTN2',evt.Value);

% === Info panel (textarea) — vyplní zbytek panelu ===
ta = uitextarea(pg,'Editable','off','FontName','Consolas','FontSize',11);
ta.Layout.Row = 13; ta.Layout.Column = [1 2];

% Ulož handlery (ať můžeš synchronizovat popisky hodnot apod.)
S.ax = ax; S.ta = ta;
S.slF = slF;   S.slN1 = slN1; S.slN2 = slN2; S.slA = slA; S.slAlpha2 = slAlpha2;
S.lblF = lblF; S.lblN1 = lblN1; S.lblN2 = lblN2; S.lblA = lblA;
S.efAlpha2 = efAlpha2; S.cbF = cbF; S.cbR = cbR; S.cbAttn2 = cbAttn2;

setappdata(uif,'State',S);






% Předkresli rozhraní a křivku
S.hIFace = plot(ax,[0 0],[-1 1],'k--','LineWidth',1); % bude se přepisovat podle YLim
S.hWave  = plot(ax, [0 0], [0 0], 'b','LineWidth',2); % placeholder
S.ax.YLimMode = 'manual';

% Popisky n1/n2 u rozhraní
S.txtN1 = text(ax,-0.02,0,'','HorizontalAlignment','right','FontSize',10,'Color',[0.2 0.2 0.2]);
S.txtN2 = text(ax, +0.02,0,'','HorizontalAlignment','left', 'FontSize',10,'Color',[0.2 0.2 0.2]);

% Callbacky (ValueChanging = plynulé přetahování)
S.slF.ValueChangingFcn  = @(src,evt) onChange('f',evt.Value);
S.slN1.ValueChangingFcn = @(src,evt) onChange('n1',evt.Value);
S.slN2.ValueChangingFcn = @(src,evt) onChange('n2',evt.Value);
S.slA.ValueChangingFcn  = @(src,evt) onChange('A0',evt.Value);
S.cbF.ValueChangedFcn   = @(src,evt) onChange('F',evt.Value);
S.cbR.ValueChangedFcn   = @(src,evt) onChange('R',evt.Value);
S.cbAttn2.ValueChangedFcn = @(src,evt) onChange('ATTN2', evt.Value);

% Sdílej stav
setappdata(uif,'State',S);

%% --- Animace (běží dokud nezavřeš okno) ---
t = 0;
while isvalid(uif)
    S = getappdata(uif,'State');

    % Odvozené veličiny
    lambda0 = S.c / S.f;  omega = 2*pi*S.f;  k0 = 2*pi / lambda0;
    v1 = S.c / S.n1;      v2 = S.c / S.n2;
    lambda1 = lambda0 / S.n1;   lambda2 = lambda0 / S.n2;
    k1 = S.n1 * k0;       k2 = S.n2 * k0;

    % Fresnel (kolmý dopad, E-pole, nemagnetická média)
    r = (S.n1 - S.n2) / (S.n1 + S.n2);
    tA = 2*S.n1 / (S.n1 + S.n2);
    R = abs(r)^2;
    T = (S.n2/S.n1) * abs(tA)^2;

    

    % Fresnel pro kolmý dopad (nemagnetická média)
    if S.showFresnell
        r  = (S.n1 - S.n2) / (S.n1 + S.n2);
        tA = 2*S.n1 / (S.n1 + S.n2);
        R  = abs(r)^2;
        T  = (S.n2/S.n1) * abs(tA)^2;
    else
        r  = 0;          % bez rozhraní žádný odraz
        tA = 1;          % plný přenos amplitudy
        R  = 0; T = 1;   % energetická bilance
    end

    % Amplitudy složek
    A_inc = S.A0;
    A_ref = (S.showReflection) * (r  * S.A0);    
    A_trn = tA * S.A0;

    % Dynamické prostorové okno (3 λ vlevo i vpravo)
    zL = linspace(-15*lambda1, 0, 600);
    zR = linspace(0,  15*lambda2, 600);
    z  = [zL, zR(2:end)];

    % Pole
    yL = A_inc * sin(omega*t - k1*zL + S.phi) + A_ref * sin(omega*t + k1*zL + S.phi);


    % Útlum v n2 (Beer–Lambert pro amplitudu)
    %if S.useAttn2
    %    yR = yR .* exp(-S.alpha2 * zR);  % amplituda ~ e^{-alpha z}
    %end

    % Přenesená vlna v n2:
    yR = A_trn * sin(omega*t - k2*zR + S.phi);
    
    % Útlum v n2 (amplituda):
    if S.useAttn2
        yR = yR .* exp(-S.alpha2 * zR);
    end

    y  = [yL, yR(2:end)];

    % Osy a křivka
    yMax = 1.2 * max([abs(A_inc)+abs(A_ref), abs(A_trn), 0.5]); % bezpečná rezerva
    S.ax.XLim = [min(z) max(z)];
    S.ax.XLim = [-1 1];
    S.ax.YLim = [-yMax yMax];
    set(S.hWave,'XData',z,'YData',y);
    set(S.hIFace,'XData',[0 0],'YData',[-yMax yMax]);
    set(S.txtN1,'Position',[-0.02*(S.ax.XLim(2)-S.ax.XLim(1)), 0], 'String',sprintf('n_1=%.3f',S.n1));
    set(S.txtN2,'Position',[+0.02*(S.ax.XLim(2)-S.ax.XLim(1)), 0], 'String',sprintf('n_2=%.3f',S.n2));

    % Infopanel
    info = sprintf([ ...
        'GLOBÁLNÍ PARAMETRY\n' ...
        '-------------------\n' ...
        'c         = %.6g\n' ...
        'f         = %.6g   (λ_0 = %.6g)\n' ...
        '\nPROSTŘEDÍ 1 (z<0)\n' ...
        '-------------------\n' ...
        'n_1       = %.6g\n' ...
        'v_1       = %.6g\n' ...
        'λ_1       = %.6g\n' ...
        'A_{inc}   = %.6g\n' ...
        'A_{ref}   = %.6g   (r = %.6g)\n' ...
        '\nPROSTŘEDÍ 2 (z≥0)\n' ...
        '-------------------\n' ...
        'n_2       = %.6g\n' ...
        'v_2       = %.6g\n' ...
        'λ_2       = %.6g\n' ...
        'A_{trn}   = %.6g   (t = %.6g)\n' ...
        '\nENERGIE\n' ...
        '-------------------\n' ...
        'R = |r|^2          = %.6g\n' ...
        'T = (n_2/n_1)|t|^2 = %.6g\n' ...
        'R + T              = %.6g\n' ...
        '\nSTAV\n' ...
        '-------------------\n' ...
        't         = %.2f s\n' ...
        ], ...
        S.c, S.f, lambda0, ...
        S.n1, v1, lambda1, S.A0, A_ref, r, ...
        S.n2, v2, lambda2, A_trn, tA, ...
        R, T, R+T, ...
        t);
    S.ta.Value = splitlines(info);

    % Aktualizuj numerické labely u sliderů (pohodlné čtení)
    S.lblF.Text  = sprintf('%.2f',S.f);
    S.lblN1.Text = sprintf('%.3f',S.n1);
    S.lblN2.Text = sprintf('%.3f',S.n2);
    S.lblA.Text  = sprintf('%.2f',S.A0);

    title(S.ax, sprintf('Šíření přes rozhraní — t = %.2f s', t));

    setappdata(uif,'State',S); % ulož zpět
    drawnow limitrate;
    pause(S.Tpause/1000);
    t = t + S.dt;
end

%% --- Lokální pomocná funkce (callbacky posuvníků/checkboxu) ---
function onChange(what, val)
    uif = gcbf; S = getappdata(uif,'State');
    switch what
        case 'f',   S.f  = max(0.10, min(5.00, val));
        case 'n1',  S.n1 = max(1.00,  min(2.50, val));
        case 'n2',  S.n2 = max(1.00,  min(2.50, val));
        case 'A0',  S.A0 = max(0.10,  min(2.00, val));
        case 'R',   S.showReflection = logical(val);
        case 'F',   S.showFresnell = logical(val);
        case 'ATTN2', S.useAttn2     = logical(val);   % << nový přepínač
        case 'ALPHA2'
            S.alpha2 = max(0, min(2.0, double(val)));   % clamp + jistota double
            % držet UI prvky synchronizované:
            if isfield(S,'slAlpha2') && isvalid(S.slAlpha2), S.slAlpha2.Value = S.alpha2; end
            if isfield(S,'efAlpha2') && isvalid(S.efAlpha2), S.efAlpha2.Value = S.alpha2; end
    end
    setappdata(uif,'State',S);
end

function out = ternary(cond,a,b)
    if cond
        out=a; 
    else
        out=b;
    end

end
