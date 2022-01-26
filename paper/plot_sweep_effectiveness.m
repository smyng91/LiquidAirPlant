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
load('paper\solution_sweep_effectiveness.mat');
Np = 20;
solution(1:2,:) = [];
mdot = solution(:,1);
eta = solution(:,2);
solution(:,23:26) = rescale(solution(:,23:26));


%%
% Pr vs T
T_hp_in = solution(:,14);
T_hp_out = solution(:,15);
T_lp_out = solution(:,16);
T_pre_Tc_out = solution(:,17);
T_rec_Tc_out = solution(:,18);
T_ptc_in = solution(:,19);
T_ptc_out = solution(:,20);
T_lp_in = solution(:,21);
T_rec_Th_out = solution(:,22);
e_pre = solution(:,23);
e_rec = solution(:,24);
e_shx1 = solution(:,25);
e_shx2 = solution(:,26);

method = 'linear';
xq = linspace(min(mdot),max(mdot),10);
x = mdot;

figure('units','normalized','position',[0.2 0.2 .5 .5])
% y = T_hp_in;
% vq = interp1(x, y, xq, method);
% h3 = plot(xq, vq, '-vk','MarkerFaceColor','k');
% hold on
% y = T_hp_out;
% vq = interp1(x, y, xq, method);
% h4 = plot(xq, vq, '-dk','MarkerFaceColor','w');
% hold on
% y = T_lp_in;
% vq = interp1(x, y, xq, method);
% h5 = plot(xq, vq, '-sk','MarkerFaceColor','w');
% hold on
y = T_lp_out;
vq = interp1(x, y, xq, method);
h6 = plot(xq, vq, '-ok','MarkerFaceColor','k');
hold on
y = T_pre_Tc_out;
vq = interp1(x, y, xq, method);
h1 = plot(xq, vq, '-+k','MarkerFaceColor','k');
hold on
y = T_rec_Tc_out;
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-^k','MarkerFaceColor','k');
hold on
y = T_rec_Th_out;
vq = interp1(x, y, xq, method);
h7 = plot(xq, vq, '-xk','MarkerFaceColor','k');
hold on
legend([h1,h2,h6,h7],'$T_3$','$T_4$','$T_8$','$T_9$','orientation','horizontal','interpreter','latex','FontSize',24,'box','off','location','northeast')

plot([0.8 0.8],[0 1000],'--r')

pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
xlim([0 2.5])
xlabel('$\dot{m}_{la}$')
ylim([100 400])
ylab=ylabel('$T~\mathrm{(K)}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.12, 0.5]);
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




% yyaxis right
% y = e_pre;
% vq = interp1(x, y, xq, method);
% h1 = plot(xq, vq, '-^k','MarkerSize',5,'MarkerFaceColor','k');
% hold on
% y = e_rec;
% vq = interp1(x, y, xq, method);
% h2 = plot(xq, vq, '-^k','MarkerSize',5,'MarkerFaceColor','k');
% hold on
% y = e_shx1;
% vq = interp1(x, y, xq, method);
% h3 = plot(xq, vq, '-^k','MarkerSize',5,'MarkerFaceColor','k');
% hold on
% y = e_shx2;
% vq = interp1(x, y, xq, method);
% h4 = plot(xq, vq, '-^k','MarkerSize',5,'MarkerFaceColor','k');
% hold on
% % ylim([0.3 0.4])
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';
% ylab=ylabel('$\eta_{ex}$','interpreter','latex');
% set(ylab, 'Units', 'Normalized', 'Position', [1.15, 0.55]);




