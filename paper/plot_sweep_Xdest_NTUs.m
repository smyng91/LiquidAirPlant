clc; clear; close all; clear fcn;
alw = 1.2;    % AxesLineWidth
fsz = 20;      % Fontsize
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
load('paper\solution_sweep_NTU.mat');
Np = 10;
solution(1:2,:) = [];
NTU_preheat = solution(:,1);
shx1 = solution(:,2);
eta = solution(:,4);
X_dest.p = solution(:,6);
X_dest.p_shx = solution(:,7);
X_dest.pre = solution(:,8);
X_dest.rec = solution(:,9);
X_dest.shx1= solution(:,10);
X_dest.shx2= solution(:,11);
X_dest.t_hp= solution(:,12);
X_dest.t_lp = solution(:,13);
X_dest.ptc = solution(:,14);
X_sun = 4.7e5;
X_in = 9.107e5;
X_u = X_sun+X_in - sum(solution(:,6:14),2);

%%
% NTU vs eta_ex
method = 'linear';
xq = linspace(min(shx1),max(shx1),20);

figure('units','normalized','position',[0.2 0.2 .5 .5])
preheat = 0.01;
x = NTU_preheat; 
y = [ X_dest.p X_dest.p_shx X_dest.pre ...
    X_dest.rec  X_dest.shx1  X_dest.shx2 ... 
    X_dest.t_hp  X_dest.t_lp  X_dest.ptc ...
    ]/(X_in+X_sun);
b=bar(x,y,'stacked','FaceColor','flat');
colormap parula
for k = 1:size(y,2)
    b(k).CData = k;
end
hold on
leg = legend(b,'$~\widetilde{X}_{dest,p_{la}}$','$~\widetilde{X}_{dest,p_{shx}}$','$~\widetilde{X}_{dest,pre}$',...
    '$~\widetilde{X}_{dest,rec}$','$~\widetilde{X}_{dest,shx,1}$','$~\widetilde{X}_{dest,shx,2}$',...
    '$~\widetilde{X}_{dest,t,1}$', '$~\widetilde{X}_{dest,t,2}$','$~\widetilde{X}_{dest,ptc}$');
set(leg,'Interpreter','latex','orientation','vertical','location','southeast','box','off')
pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
xlim([0 9])
xlabel('$\mathrm{NTU}_{pre}$')
ylim([0 .4])
ylab=ylabel('$\dot{X}_{dest,i}/\dot{X}_{in}$');
set(ylab, 'Units', 'Normalized', 'Position', [-0.12, 0.55]);
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
% title('$\dot{m}_{la}=0.5~\mathrm{kg/s}$')
yyaxis right
y = eta;
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-^k','MarkerSize',5,'MarkerFaceColor','k');
ylim([0.45 0.55])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ylab=ylabel('$\eta_{ex}$','interpreter','latex');
set(ylab, 'Units', 'Normalized', 'Position', [1.15, 0.55]);


%%
% NTU vs T
T_hp_in = solution(:,15);
T_hp_out = solution(:,16);
T_lp_out = solution(:,17);
T_pre_Tc_out = solution(:,18);
T_rec_Tc_out = solution(:,19);
T_ptc_in = solution(:,20);
T_ptc_out = solution(:,21);
T_lp_in = solution(:,22);

method = 'linear';
% xq = linspace(min(p_r),max(p_r),20);

figure('units','normalized','position',[0.2 0.2 .5 .5])
mdot_la = 0.5;
y = T_hp_in;
vq = interp1(x, y, xq, method);
h3 = plot(xq, vq, '-vk','MarkerFaceColor','k');
hold on
 y = T_hp_out;
vq = interp1(x, y, xq, method);
h4 = plot(xq, vq, '-dk','MarkerFaceColor','w');
hold on
 y = T_lp_in;
vq = interp1(x, y, xq, method);
h5 = plot(xq, vq, '-sk','MarkerFaceColor','w');
hold on
 y = T_lp_out;
vq = interp1(x, y, xq, method);
h6 = plot(xq, vq, '-ok','MarkerFaceColor','w');
hold on
 y = T_pre_Tc_out;
vq = interp1(x, y, xq, method);
h1 = plot(xq, vq, '-+k','MarkerFaceColor','k');
hold on
 y = T_rec_Tc_out;
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-^k','MarkerFaceColor','k');
hold on
legend([h1,h2,h3,h4,h5,h6],'$T_3$','$T_4$','$T_5$','$T_6$','$T_7$','$T_8$','orientation','horizontal','interpreter','latex','FontSize',24,'box','off','location','south')

pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
% xlim([-1 6])
xlabel('$\mathrm{NTU}_{pre}$')
ylim([0 700])
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

