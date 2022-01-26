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
load('paper\solution_sweep_eta_rt.mat');
Np = 10;
mdot = solution(:,1);
p_r = solution(:,2);
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
% Pr vs eta_ex
method = 'pchip';
xq = linspace(min(p_r),max(p_r),20);

figure('units','normalized','position',[0.2 0.2 .5 .5])
mdot_la = 0.5;
x = p_r(abs(mdot-mdot_la)<1e-1); 
y = [ X_dest.p(abs(mdot-mdot_la)<1e-1) X_dest.p_shx(abs(mdot-mdot_la)<1e-1) X_dest.pre(abs(mdot-mdot_la)<1e-1) ...
    X_dest.rec(abs(mdot-mdot_la)<1e-1)  X_dest.shx1(abs(mdot-mdot_la)<1e-1)  X_dest.shx2(abs(mdot-mdot_la)<1e-1) ... 
    X_dest.t_hp(abs(mdot-mdot_la)<1e-1)  X_dest.t_lp(abs(mdot-mdot_la)<1e-1)  X_dest.ptc(abs(mdot-mdot_la)<1e-1) ...
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
set(leg,'Interpreter','latex','orientation','vertical','location','northwest','box','off')
pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
xlim([-5 60])
xlabel('$\tilde{p}_r$')
ylim([0 1])
ylab=ylabel('$\dot{X}_{dest,i}/\dot{X}_{in}$');
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
title('$\dot{m}_{la}=0.5~\mathrm{kg/s}$')
yyaxis right
y = eta(abs(mdot-mdot_la)<1e-1);
vq = interp1(x, y, xq, method);
h2 = plot(xq, vq, '-^k','MarkerSize',5,'MarkerFaceColor','k');
% ylim([0.3 0.4])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ylab=ylabel('$\eta_{ex}$','interpreter','latex');
set(ylab, 'Units', 'Normalized', 'Position', [1.15, 0.55]);


