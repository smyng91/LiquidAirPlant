clc; clear; close all; clear fcn;
alw = 1.2;    % AxesLineWidth
fsz = 18;      % Fontsize
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
set(0,'defaulttextinterpreter','latex');
mrk_size = 3;

%%
clc;
load('paper\solution_sensitivity.mat');
xlabel = {'$\dot{m}_{la}$';'$\dot{m}_{aa}$';'$\dot{m}_{ptc}$';'$\dot{m}_{shx,1}$';'$\tilde{p}_r$';...
    '$\mathrm{NTU}_{pre}$';'$\mathrm{NTU}_{rec}$';'$\mathrm{NTU}_{shx,1}$';'$\mathrm{NTU}_{shx,2}$'};
n = 25;
p2 = 20e6;
p0 = 0.1e6;
p_int = p0 + (p2-p0)* solution(:,3);
p_r = p2./p_int;
mdot_shx1 = solution(:,4).*solution(:,5);
x = [ solution(:,1), solution(:,2),solution(:,4), mdot_shx1, p_r, solution(:,6), solution(:,7), solution(:,8), solution(:,9) ];

%%
nt = 10; % 10 - eta_rt; 11 - eta_II; 12 - COP
ymin = linspace(min(solution(:,nt)),max(solution(:,nt)),n);
group = ones(length(solution(:,nt)),1);
for i=1:1:n - 1
    group(solution(:,nt)>=ymin(i) & solution(:,nt)<=ymin(i+1)) = i+1;
end
colorMap = jet(n);
figure('units','normalized','position',[0 0 .9 .9])
gplotmatrix(x,[],group,[colorMap],['x'],mrk_size,'off','grpbars',xlabel,xlabel)
hold on
% remove upper triangle
% k = 0;
% for i = 1:9
%     for j = i+1:9
%         delete(subplot(9,9,j+k));
%     end
%     k = k+9;
% end
pbaspect([1 1 1])
sgtitle('Round trip efficiency ($\eta_{rt}$)','FontSize',30)
colormap(colorMap)
caxis([ymin(1) ymin(end)])
clb=colorbar;
clb.YTick = ymin;
clb.YTickLabel = (sprintf('%2.3g\n',linspace(ymin(1),ymin(end),length(solution(1,:)))));
colorbar('Position', [.91 0.11  0.01  0.81])
print('paper\mplot_eta_rt','-djpeg','-r400')



%%
nt = 11; % 10 - eta_rt; 11 - eta_II; 12 - COP
ymin = linspace(min(solution(:,nt)),max(solution(:,nt)),n);
group = ones(length(solution(:,nt)),1);
for i=1:1:n - 1
    group(solution(:,nt)>=ymin(i) & solution(:,nt)<=ymin(i+1)) = i+1;
end
colorMap = jet(n);
figure('units','normalized','position',[0 0 .9 .9])
gplotmatrix(x,[],group,[colorMap],['x'],mrk_size,'off','grpbars',xlabel,xlabel)
hold on
% remove upper triangle
% k = 0;
% for i = 1:9
%     for j = i+1:9
%         delete(subplot(9,9,j+k));
%     end
%     k = k+9;
% end
pbaspect([1 1 1])
sgtitle('Exergy efficiency ($\eta_{ex}$)','FontSize',30)
colormap(colorMap)
caxis([ymin(1) ymin(end)])
clb=colorbar;
clb.YTick = ymin;
clb.YTickLabel = (sprintf('%2.3g\n',linspace(ymin(1),ymin(end),length(solution(1,:)))));
colorbar('Position', [.91 0.11  0.01  0.81])
print('paper\mplot_eta_ex','-djpeg','-r400')


