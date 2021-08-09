clc; clear; close all; clear fcn;
alw = 1.2;    % AxesLineWidth
fsz = 24;      % Fontsize
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

%%
clc;
load('paper\solution_sweep_eta_rt.mat');
Np = 10;
mdot = solution(:,1);
p_r = solution(:,2);
eta = solution(:,3);

%%
% Pr vs eta_rt
method = 'pchip';
xq = linspace(min(p_r),max(p_r),100);

figure('units','normalized','position',[0.2 0.2 .5 .5])
mdot_la = 0.5;
x = p_r(abs(mdot-mdot_la)<1e-1); y = eta(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h1 = plot(xq, vq, '-k');
hold on
mdot_la = 0.9;
x = p_r(abs(mdot-mdot_la)<1e-1); y = eta(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-k');
mdot_la = 1.5;
x = p_r(abs(mdot-mdot_la)<1e-1); y = eta(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h3 = plot(xq, vq, '-k');
mdot_la = 1.7;
x = p_r(abs(mdot-mdot_la)<1e-1); y = eta(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h4 = plot(xq, vq, '-k');

pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
% xlim([-1 6])
xlabel('$\tilde{p}_r$')
% ylim([0 0.25])
ylab=ylabel('$\eta_{rt}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.55]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',2*get(gca,'ticklength'))
hold on

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
set(b,'LineWidth',alw)
pbaspect([1.2 1 1])
axes(a)


%%
% Pr vs T
T_hp_in = solution(:,15);
T_hp_out = solution(:,16);
T_lp_out = solution(:,17);
T_pre_Tc_out = solution(:,18);
T_rec_Tc_out = solution(:,19);
T_ptc_in = solution(:,20);
T_ptc_out = solution(:,21);
T_lp_in = solution(:,22);

method = 'pchip';
xq = linspace(min(p_r),max(p_r),20);

figure('units','normalized','position',[0.2 0.2 .5 .5])
mdot_la = 0.5;
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_hp_in(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h3 = plot(xq, vq, '-vk','MarkerFaceColor','k');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_hp_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h4 = plot(xq, vq, '-dk','MarkerFaceColor','w');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_lp_in(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h5 = plot(xq, vq, '-sk','MarkerFaceColor','w');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_lp_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h6 = plot(xq, vq, '-ok','MarkerFaceColor','w');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_pre_Tc_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h1 = plot(xq, vq, '-+k','MarkerFaceColor','k');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_rec_Tc_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-^k','MarkerFaceColor','k');
hold on
legend([h1,h2,h3,h4,h5,h6],'$T_3$','$T_4$','$T_5$','$T_6$','$T_7$','$T_8$','orientation','horizontal','interpreter','latex','FontSize',24,'box','off','location','south')

pbaspect([1.2 1 1])
title('$\dot{m}_{la}=0.5~\mathrm{kg/s}$','interpreter','latex')
set(gca,'LineWidth',alw)
% xlim([-1 6])
xlabel('$\tilde{p}_{r}$')
ylim([200 700])
ylab=ylabel('$T~\mathrm{(K)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',2*get(gca,'ticklength'))
hold on

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
set(b,'LineWidth',alw)
pbaspect([1.2 1 1])
axes(a)


figure('units','normalized','position',[0.2 0.2 .5 .5])
mdot_la = 1.5;
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_hp_in(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h3 = plot(xq, vq, '-vk','MarkerFaceColor','k');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_hp_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h4 = plot(xq, vq, '-dk','MarkerFaceColor','w');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_lp_in(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h5 = plot(xq, vq, '-sk','MarkerFaceColor','w');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_lp_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h6 = plot(xq, vq, '-ok','MarkerFaceColor','w');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_pre_Tc_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h1 = plot(xq, vq, '-+k','MarkerFaceColor','k');
hold on
x = p_r(abs(mdot-mdot_la)<1e-1); y = T_rec_Tc_out(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-^k','MarkerFaceColor','k');
hold on
legend([h1,h2,h3,h4,h5,h6],'$T_3$','$T_4$','$T_5$','$T_6$','$T_7$','$T_8$','orientation','horizontal','interpreter','latex','FontSize',24,'box','off','location','south')

pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
% xlim([-1 6])
xlabel('$\tilde{p}_{r}$')
ylim([0 500])
ylab=ylabel('$T~\mathrm{(K)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',2*get(gca,'ticklength'))
title('$\dot{m}_{la}=1.5~\mathrm{kg/s}$','interpreter','latex')
hold on

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
set(b,'LineWidth',alw)
pbaspect([1.2 1 1])
axes(a)

