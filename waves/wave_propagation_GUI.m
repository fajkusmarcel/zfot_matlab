function wave_propagation_GUI()
% Jednoduché GUI: postupná harmonická vlna v JEDNOM prostředí
% Ovládání: amplituda A, frekvence f, index lomu n, délka osy z_max
% Info panel: c, n, v, f, lambda, k, omega, A, t
%
% Pozn.: používáme normalizovanou rychlost c = 1 (libovolné jednotky).

%% --------- Stav (výchozí hodnoty) ----------
S.c      = 1.0;     % "rychlost světla ve vakuu" v norm. jednotkách
S.f      = 1.0;     % frekvence [Hz]
S.A      = 1.0;     % amplituda
S.n      = 1.5;     % index lomu prostředí
S.zmax   = 10.0;    % délka osy z [0..zmax]
S.dt     = 0.02;    % časový krok animace [s]
S.tpause = 0.02;    % pauza mezi snímky [s]
S.t      = 0.0;     % čas
S.run    = false;   % běží animace?

%% --------- Okno a layout ----------
fig = figure('Name','Vlna v jednom prostredi','Color','w',...
    'Units','pixels','Position',[100 100 1000 600]);
set(fig,'CloseRequestFcn',@onClose);

% Levé osy (graf)
ax = axes('Parent',fig,'Units','pixels','Position',[60 80 600 480]);
hold(ax,'on'); box(ax,'on'); grid(ax,'on'); ax.GridAlpha = 0.25;
xlabel(ax,'z'); ylabel(ax,'E (arb. u.)'); title(ax,'Postupna harmonicka vlna (1 prostredi)');

% Pravý panel (ovladače + info)
p = uipanel('Parent',fig,'Title','Parametry a informace','Units','pixels',...
    'Position',[700 60 280 510],'BackgroundColor','w');

%% --------- Křivka vlny ----------
z = linspace(0, S.zmax, 1200);
hWave = plot(ax, z, zeros(size(z)),'LineWidth',2);

%% --------- Ovládací prvky ----------
y0 = 460; dy = 60; valW = 60; sW = 180; x1 = 10; x2 = 110;

% Frekvence
uicontrol(p,'Style','text','String','f [Hz]','HorizontalAlignment','left',...
    'Units','pixels','Position',[x1 y0 80 20],'BackgroundColor','w','FontWeight','bold');
txtF = uicontrol(p,'Style','text','String',num2str(S.f,'%.2f'),...
    'Units','pixels','Position',[x2 y0 valW 20],'BackgroundColor','w','HorizontalAlignment','left');
uicontrol(p,'Style','slider','Min',1,'Max',5,'Value',S.f,...
    'Units','pixels','Position',[x1 y0-18 sW 18],'Callback',@(u,~) setF(u.Value));

% Amplituda
y = y0-dy;
uicontrol(p,'Style','text','String','A','HorizontalAlignment','left',...
    'Units','pixels','Position',[x1 y 80 20],'BackgroundColor','w','FontWeight','bold');
txtA = uicontrol(p,'Style','text','String',num2str(S.A,'%.2f'),...
    'Units','pixels','Position',[x2 y valW 20],'BackgroundColor','w','HorizontalAlignment','left');
uicontrol(p,'Style','slider','Min',0.1,'Max',1,'Value',S.A,...
    'Units','pixels','Position',[x1 y-18 sW 18],'Callback',@(u,~) setA(u.Value));

% Index lomu
y = y - dy;
uicontrol(p,'Style','text','String','n (index lomu)','HorizontalAlignment','left',...
    'Units','pixels','Position',[x1 y 140 20],'BackgroundColor','w','FontWeight','bold');
txtN = uicontrol(p,'Style','text','String',num2str(S.n,'%.3f'),...
    'Units','pixels','Position',[x2 y valW 20],'BackgroundColor','w','HorizontalAlignment','left');
uicontrol(p,'Style','slider','Min',1.0,'Max',2.5,'Value',S.n,...
    'Units','pixels','Position',[x1 y-18 sW 18],'Callback',@(u,~) setN(u.Value));

% Z-osa max (0..X)
y = y - dy;
uicontrol(p,'Style','text','String','z_{max}','HorizontalAlignment','left',...
    'Units','pixels','Position',[x1 y 80 20],'BackgroundColor','w','FontWeight','bold');
txtZ = uicontrol(p,'Style','text','String',num2str(S.zmax,'%.2f'),...
    'Units','pixels','Position',[x2 y valW 20],'BackgroundColor','w','HorizontalAlignment','left');
uicontrol(p,'Style','slider','Min',1,'Max',10,'Value',S.zmax,...
    'Units','pixels','Position',[x1 y-18 sW 18],'Callback',@(u,~) setZ(u.Value));

% Start/Stop
y = y - dy + 10;
btn = uicontrol(p,'Style','togglebutton','String','Start','Units','pixels',...
    'Position',[x1 y 100 28],'Callback',@onRun);

% Info pole (READ-ONLY)
yInfo = 20;
info = uicontrol(p,'Style','edit','Max',2,'Min',0,'Enable','inactive',...
    'HorizontalAlignment','left','BackgroundColor',[0.97 0.97 0.97],...
    'Units','pixels','Position',[10 yInfo 260 180],'FontName','Consolas');

%% --------- Pomocné (update funkcionalita) ----------
    function recompute_and_redraw()
        % Osové body
        z = linspace(0, S.zmax, 1200);
        % Odvozené veličiny
        v  = S.c / S.n;                         % fazova rychlost v prostredi
        if S.f <= 0
            lambda = Inf;
        else
            lambda = v / S.f;
        end
        if isfinite(lambda) && lambda>0
            k  = 2*pi / lambda;                 % vlnove cislo
        else
            k  = 0;
        end
        omega = 2*pi*S.f;                       % uhlova frekvence

        % Vlna (postupna)
        yWave = S.A * sin(omega*S.t - k*z);

        % Osy a kresba
        ax.XLim = [0 S.zmax];
        yAbs = max(1.2*max(1e-6, abs(S.A)), 0.5);
        ax.YLim = [-1 1];
        set(hWave,'XData',z,'YData',yWave);

        % Info text (čisté ASCII, žádné uvozovky "")
        txt = sprintf([ ...
            'ZAKLADNI PARAMETRY\n' ...
            '-------------------\n' ...
            'c          = %.6g\n' ...
            'n          = %.6g\n' ...
            'v=c/n      = %.6g\n' ...
            'f          = %.6g\n' ...
            'lambda=v/f = %.6g\n' ...
            'k=2*pi/l   = %.6g\n' ...
            'omega=2*pi*f = %.6g\n' ...
            'A          = %.6g\n' ...
            'z_max      = %.6g\n' ...
            '\n' ...
            't          = %.2f s\n' ], ...
            S.c, S.n, v, S.f, lambda, k, omega, S.A, S.zmax, S.t);
        set(info,'String',txt);

        % Titulek
        title(ax, sprintf('Postupna vlna v homogenim prostredi — t = %.2f s', S.t));
        drawnow limitrate;
    end

    function setF(val)
        S.f = max(0.05, min(5.0, val));
        set(txtF,'String',num2str(S.f,'%.2f'));
        recompute_and_redraw();
    end

    function setA(val)
        S.A = max(0.1, min(2.0, val));
        set(txtA,'String',num2str(S.A,'%.2f'));
        recompute_and_redraw();
    end

    function setN(val)
        S.n = max(1.0, min(2.5, val));
        set(txtN,'String',num2str(S.n,'%.3f'));
        recompute_and_redraw();
    end

    function setZ(val)
        S.zmax = max(1, min(50, val));
        set(txtZ,'String',num2str(S.zmax,'%.2f'));
        recompute_and_redraw();
    end

    function onRun(src,~)
        S.run = logical(get(src,'Value'));
        set(src,'String', tern(S.run,'Stop','Start'));
        while ishghandle(fig) && S.run
            S.t = S.t + S.dt;
            recompute_and_redraw();
            pause(S.tpause);
            S.run = logical(get(src,'Value')); % respektuj kliknuti behem behu
        end
        if ishghandle(src)
            set(src,'String', tern(S.run,'Stop','Start'));
        end
    end

    function onClose(~,~)
        S.run = false;
        if ishghandle(fig)
            delete(fig);
        end
    end

    function out = tern(cond,a,b), if cond, out=a; else, out=b; end, end

%% --------- První vykreslení ----------
recompute_and_redraw();

end
