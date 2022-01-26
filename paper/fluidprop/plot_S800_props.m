clc; clear; close all; clear fcn;
alw = 1.5;    % AxesLineWidth
fsz = 16;      % Fontsize
ftype = 'Times New Roman';  % Font type
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
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

T = linspace(300,650,100);
p = 0.1e6;

for i=1:length(T)
        k(i) = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','P',p,'T',T(i),'INCOMP::S800');
        cp(i) = py.CoolProp.CoolProp.PropsSI('CPMASS','P',p,'T',T(i),'INCOMP::S800');
        mu(i) = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',p,'T',T(i),'INCOMP::S800');
        rho(i) = py.CoolProp.CoolProp.PropsSI('Dmass','P',p,'T',T(i),'INCOMP::S800');
end

figure('units','normalized','position',[0.2 0.2 .2 .5])
sgtitle('Syltherm 800','FontSize',20)
subplot(2,2,1)
scatter(T(1:10:end),k(1:10:end),40,'ok','filled')
ylab=ylabel('$k~\mathrm{(W/m\cdot K)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
ylim([0.05, 0.15])
xlabel('$T~\mathrm{(K)}$')
% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
axes(a)
hold on
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.15 * (yl(2)-yl(1)) + yl(1);
p = polyfit(T,k,1);
caption = sprintf('$k = %f T +$\n $~~~~%f$', p(1), p(2));
text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
f1 = polyval(p,T);
plot(T,f1,'r-')

subplot(2,2,2)
scatter(T(1:10:end),cp(1:10:end)/1e3,40,'ok','filled')
xlabel('$T~\mathrm{(K)}$')
ylab=ylabel('$c_p~\mathrm{(kJ/kg\cdot K)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
axes(a)
hold on
xl = xlim;
yl = ylim;
xt = 0.55 * (xl(2)-xl(1)) + xl(1);
yt = 0.15 * (yl(2)-yl(1)) + yl(1);
p = polyfit(T,cp/1e3,1);
caption = sprintf('$c_p = %f T +$\n $~~~~%f$', p(1), p(2));
text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
f1 = polyval(p,T);
plot(T,f1,'r-')


subplot(2,2,3)
scatter(T(1:10:end),mu(1:10:end)*1e3,40,'ok','filled')
xlabel('$T~\mathrm{(K)}$')
ylab=ylabel('$\mu~\mathrm{(mPa\cdot s)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
axes(a)
hold on
xl = xlim;
yl = ylim;
xt = 0.2 * (xl(2)-xl(1)) + xl(1);
yt = 0.85 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('$\\mu/1000 = 17.15e^{-0.027T}+$ \n$~~~~~~~~~~~~0.041e^{-0.0079T}$');
text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
f1 = @(x) 17.15*exp(-0.02671*x) + 0.04084*exp(-0.007925*x);
plot(T,f1(T)*1e3,'r-')


subplot(2,2,4)
scatter(T(1:10:end),rho(1:10:end),40,'ok','filled')
xlabel('$T~\mathrm{(K)}$')
ylab=ylabel('$\rho~\mathrm{(kg/m^3)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
axes(a)
hold on
xl = xlim;
yl = ylim;
xt = 0.45 * (xl(2)-xl(1)) + xl(1);
yt = 0.75 * (yl(2)-yl(1)) + yl(1);
p = polyfit(T,rho,2);
caption = sprintf(['$\\rho = %f T^2 $\n $~~~~~~%f T$ +' ...
    '\n$~~~~~~~~~%f$'], p(1), p(2), p(3));
text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
f1 = polyval(p,T);
plot(T,f1,'r-')
