alw = 1.5;    % AxesLineWidth
fsz = 25;      % Fontsize
ftype = 'Times New Roman';  % Font type
lw = 2;      % LineWidth
msz = 100;       % MarkerSize
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultUicontrolFontName',ftype);
set(0,'defaultUitableFontName',ftype);
set(0,'defaultAxesFontName',ftype);
set(0,'defaultTextFontName',ftype);
set(0,'defaultUipanelFontName',ftype);
set(0,'defaultAxesFontSize',fsz)
set(0,'defaulttextfontsize',fsz);
set(0,'defaulttextinterpreter','latex')
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

 global data_Ts_AA data_Ts_AANG data_Ts_RAANG

 
%% Baseline
fluid = 'AIR';
% create saturation dome
Tcrit = py.CoolProp.CoolProp.PropsSI('Tcrit','',0,'',0,fluid);
T = linspace(0,Tcrit,1000);

for i=1:length(T)
    try
        s_liq(i) = py.CoolProp.CoolProp.PropsSI('Smass','Q',0,'T',T(i),fluid);
        s_gas(i) = py.CoolProp.CoolProp.PropsSI('Smass','Q',1,'T',T(i),fluid);
        T_sat(i) = T(i);
    catch
        s_liq(i) = NaN;
        s_gas(i) = NaN;
        T_sat(i) = NaN;
        continue
    end
end
% obtain isobars
T = linspace(0,2000,500);
p = [0.1, 0.2, 1.5, 3, 20]*1e6;
for i=1:length(p)
    for j=1:length(T)
        try
            s(j,i) = py.CoolProp.CoolProp.PropsSI('Smass','P',p(i),'T',T(j),fluid);
            T_iso(j,i) = T(j);
        catch
            s(j,i) = NaN;
            T_iso(j,i) = NaN;
            continue
        end
    end
end

%%
data = data_Ts_AA;
[n, ~] = size(data);
figure('units','normalized','position',[0.1 0.2 1 0.5])
subplot(1,3,1)
plot(s_liq/1e3,T_sat,'-k')
hold on
plot(s_gas/1e3,T_sat,'-k')
plot(s/1e3,T,'--k')

for i=[4,6]
    plot([data(i-1,1)/1e3,data(i,1)/1e3],[data(i-1,2),data(i,2)],'-b')
end

plot(s(s(:,5)>=data(2,1) & s(:,5)<=data(3,1),5)/1e3, T_iso(T_iso(:,5)>=data(2,2) & T_iso(:,5)<=data(3,2),5),'-b')
plot(s(s(:,4)>=data(4,1) & s(:,4)<=data(5,1),4)/1e3, T_iso(T_iso(:,4)>=data(4,2) & T_iso(:,4)<=data(5,2),4),'-b')
scatter(data(:,1)/1e3,data(:,2),msz,'or','MarkerFaceColor','r','MarkerEdgeColor','k')

fs=15;
h=text(2.4,400,'$20~\mathrm{MPa}$','FontSize',fs);
set(h,'Rotation',75);
h=text(3,400,'$3~\mathrm{MPa}$','FontSize',fs);
set(h,'Rotation',75);
h=text(3.3,270,'$1.5~\mathrm{MPa}$','FontSize',fs);
set(h,'Rotation',75);
h=text(3.75,380,'$0.2~\mathrm{MPa}$','FontSize',fs);
set(h,'Rotation',75);
h=text(4.35,380,'$0.1~\mathrm{MPa}$','FontSize',fs);
set(h,'Rotation',75);
h=text(-0.5,450,'AA');
px = data(1:n,1)/1e3;
py = data(1:n,2);
for i=1:length(px)
    str2 = sprintf('%i',i);
    text(px(i)-.4, py(i), str2, 'Interpreter', 'latex')
end

pbaspect([1.4 1 1])
set(gca,'LineWidth',alw)
xlim([-1 6])
xlabel('$s~(\mathrm{kJ/kg\cdot K})$')
ylim([0 500])
ylab=ylabel('$T~(\mathrm{K})$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))
hold on

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
set(b,'LineWidth',alw)
pbaspect([1.4 1 1])
axes(a)
% print('Ts_baseline','-djpeg','-r400')


%%
data = data_Ts_AANG;

subplot(1,3,2)
plot(s_liq/1e3,T_sat,'-k')
hold on
plot(s_gas/1e3,T_sat,'-k')
plot(s/1e3,T,'--k')

plot([data(1,1)/1e3,data(2,1)/1e3],[data(1,2),data(2,2)],'-b')
plot([data(5,1)/1e3,data(6,1)/1e3],[data(5,2),data(6,2)],'-b')
plot([data(7,1)/1e3,data(8,1)/1e3],[data(7,2),data(8,2)],'-b')
i=2; j=5;
plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
i=3; j=5;
plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
i=4; j=5;
plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
i=6; j=4;
plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
i=8; j=1;
plot(s(s(:,j)<=data(i,1) & s(:,j)>=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)<=data(i,2) & T_iso(:,j)>=data(i+1,2),j),'-b')

scatter(data(:,1)/1e3,data(:,2),msz,'or','MarkerFaceColor','r','MarkerEdgeColor','k')


h=text(3.5,80,'AANG');
px = data(:,1)/1e3;
py = data(:,2);
for i=1:length(px)
    str2 = sprintf('%i',i);
    text(px(i)-.4, py(i), str2, 'Interpreter', 'latex')
end

pbaspect([1.4 1 1])
set(gca,'LineWidth',alw)
xlim([-1 6])
xlabel('$s~(\mathrm{kJ/kg\cdot K})$')
ylim([0 2000])
ylab=ylabel('$T~(\mathrm{K})$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.18, 0.5]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))
hold on

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
set(b,'LineWidth',alw)
pbaspect([1.4 1 1])
axes(a)
% print('Ts_baseline','-djpeg','-r400')

 
%% Cold Brayton
data = data_Ts_RAANG;

subplot(1,3,3)
plot(s_liq/1e3,T_sat,'-k')
hold on
plot(s_gas/1e3,T_sat,'-k')
plot(s/1e3,T,'--k')

plot([data(1,1)/1e3,data(2,1)/1e3],[data(1,2),data(2,2)],'-b')
plot([data(5,1)/1e3,data(6,1)/1e3],[data(5,2),data(6,2)],'-b')
% i=2; j=5;
% plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
% i=3; j=5;
% plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
% i=4; j=5;
% plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
% i=6; j=4;
% plot(s(s(:,j)>=data(i,1) & s(:,j)<=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)>=data(i,2) & T_iso(:,j)<=data(i+1,2),j),'-b')
% i=8; j=1;
% plot(s(s(:,j)<=data(i,1) & s(:,j)>=data(i+1,1),j)/1e3, T_iso(T_iso(:,j)<=data(i,2) & T_iso(:,j)>=data(i+1,2),j),'-b')

scatter(data(:,1)/1e3,data(:,2),msz,'or','MarkerFaceColor','r','MarkerEdgeColor','k')

h=text(3.5,80,'RAANG');
px = [data(:,1);data(end,1)]/1e3;
py = [data(:,2);data(end,2)];
for i=1:length(px)
    str2 = sprintf('%ia',i);
    text(px(i)-.4, py(i), str2, 'Interpreter', 'latex')
end

pbaspect([1.4 1 1])
set(gca,'LineWidth',alw)
xlim([-1 6])
xlabel('$s~(\mathrm{kJ/kg\cdot K})$')
ylim([0 2000])
ylab=ylabel('$T~(\mathrm{K})$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.18, 0.5]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))
hold on

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
set(b,'LineWidth',alw)
pbaspect([1.4 1 1])
axes(a)
% print('Ts_baseline','-djpeg','-r400')
