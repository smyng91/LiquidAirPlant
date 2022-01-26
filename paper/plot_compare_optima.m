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
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  



eta_rt = [6.9120    0.5707  5.1406    0.5214    2.9179    2.8686    2.7397    3.4738    0.1068  ]';
eta_ex = [0.3312    0.2266    9.5124    0.4933    2.2631    1.7217    3.4825    4.5327  0.1445    ]';
balance = [...
    1.076213198091668;      
    0.548587299087344;  
    4.172148090860937;  
    0.155894182073230;  
    3.739407410183461;  
    2.249099165017271;  
    4.327127891503268;  
    1.685345533978621
    0.122806279847624;  
];

A = [eta_rt, eta_ex, balance];
B = A ./ A(:,3);

figure('units','normalized','position',[0.2 0.2 .5 .5])
barh(B, 1)
pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
% xlim([-1 6])
xlabel('Relative value')
% ylim([0 0.25])
ylab=ylabel('');
set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.55]);
set(gca,'TickDir','out','YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',2*get(gca,'ticklength'))
hold on
yticklabels({'$\dot{m}_{aa}$' '$\dot{m}_{la}$' '$\dot{m}_{ptc}$' '$\dot{m}_{shx,1}$'...
    '$\mathrm{NTU}_{pre}$' '$\mathrm{NTU}_{rec}$' '$\mathrm{NTU}_{shx,1}$' '$\mathrm{NTU}_{shx,2}$' '$\tilde{p}_{r}$'})
legend({'$~$Maximum $\eta_{rt}$','$~$Maximum $\eta_{ex}$','$~$Balanced (Baseline)'},'box','off')

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
