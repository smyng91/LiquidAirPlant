%%
% plot T-s diagram for AAS configuration

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

fluid = PLANT.FLUID;

%%
T = [
    PUMP.T_in, PUMP.T_out;
    PREHEATER.T_c_in, PREHEATER.T_c_out;
    RECUPERATOR.T_c_in, RECUPERATOR.T_c_out;
    SHX1.T_c_in, SHX1.T_c_out;
    TURBINE_HP.T_in, TURBINE_HP.T_out;
    SHX2.T_c_in, SHX2.T_c_out;
    TURBINE_LP.T_in, TURBINE_LP.T_out;
    RECUPERATOR.T_h_in, RECUPERATOR.T_h_out
];
smass = @(T,p,fluid) py.CoolProp.CoolProp.PropsSI('Smass','T',T,'P',p,fluid);
s = [
    smass(T(1,1), PUMP.p_in, fluid), smass(T(1,2), PUMP.p_out, fluid) 
    smass(T(2,1), PREHEATER.p_c_in, fluid), smass(T(2,2), PREHEATER.p_c_in, fluid) 
    smass(T(3,1), RECUPERATOR.p_c_in, fluid), smass(T(3,2), RECUPERATOR.p_c_in, fluid) 
    smass(T(4,1), SHX1.p_c_in, fluid), smass(T(4,2), SHX1.p_c_in, fluid) 
    smass(T(5,1), TURBINE_HP.p_in, fluid), smass(T(5,2), TURBINE_HP.p_out, fluid) 
    smass(T(6,1), SHX2.p_c_in, fluid), smass(T(6,2), SHX2.p_c_in, fluid) 
    smass(T(7,1), TURBINE_LP.p_in, fluid), smass(T(7,2), TURBINE_LP.p_out, fluid) 
    smass(T(8,1), RECUPERATOR.p_h_in, fluid), smass(T(8,2), RECUPERATOR.p_h_in, fluid)    
    ];
p = [PUMP.p_in, PUMP.p_out, TURBINE_HP.p_out, TURBINE_LP.p_out];


% create saturation dome
Tcrit = py.CoolProp.CoolProp.PropsSI('Tcrit','',0,'',0,fluid);
Tsat = linspace(0,Tcrit,2000);

for i=1:length(Tsat)
    try
        s_liq(i) = py.CoolProp.CoolProp.PropsSI('Smass','Q',0,'T',Tsat(i),fluid);
        s_gas(i) = py.CoolProp.CoolProp.PropsSI('Smass','Q',1,'T',Tsat(i),fluid);
        T_sat(i) = Tsat(i);
    catch
        s_liq(i) = NaN;
        s_gas(i) = NaN;
        T_sat(i) = NaN;
        continue
    end
end
% obtain isobars
Tsat = linspace(0,2000,500);

for i=1:length(p)
    for j=1:length(Tsat)
        try
            s_iso(j,i) = py.CoolProp.CoolProp.PropsSI('Smass','P',p(i),'T',Tsat(j),fluid);
            T_iso(j,i) = Tsat(j);
        catch
            s_iso(j,i) = NaN;
            T_iso(j,i) = NaN;
            continue
        end
    end
end

figure('units','normalized','position',[0.2 0.2 0.5 0.5])
plot(s_liq/1e3,T_sat,'-k')
hold on
plot(s_gas/1e3,T_sat,'-k')
plot(s_iso/1e3,Tsat,'--k')
plot( [s(1,1)/1e3,s(1,2)/1e3],[T(1,1),T(1,2)], '-b' )
plot( s_iso( s_iso(:,2)>=s(2,1) & s_iso(:,2)<=s(2,2),2)/1e3, T_iso(T_iso(:,2)>=T(2,1) & T_iso(:,2)<=T(2,2),2),'-b')
plot( s_iso( s_iso(:,2)>=s(3,1) & s_iso(:,2)<=s(3,2),2)/1e3, T_iso(T_iso(:,2)>=T(3,1) & T_iso(:,2)<=T(3,2),2),'-b')
plot( s_iso( s_iso(:,2)>=s(4,1) & s_iso(:,2)<=s(4,2),2)/1e3, T_iso(T_iso(:,2)>=T(4,1) & T_iso(:,2)<=T(4,2),2),'-b')
plot( [s(5,1)/1e3,s(5,2)/1e3],[T(5,1),T(5,2)], '-b' )
plot( s_iso( s_iso(:,3)>=s(6,1) & s_iso(:,3)<=s(6,2),3)/1e3, T_iso(T_iso(:,3)>=T(6,1) & T_iso(:,3)<=T(6,2),3),'-b')
plot( [s(7,1)/1e3,s(7,2)/1e3],[T(7,1),T(7,2)], '-b' )
plot( s_iso( s_iso(:,4)<=s(8,1) & s_iso(:,4)>=s(8,2),4)/1e3, T_iso(T_iso(:,4)<=T(8,1) & T_iso(:,4)>=T(8,2),4),'-b')
scatter(s(:,1)/1e3,T(:,1),msz,'or','MarkerFaceColor','r','MarkerEdgeColor','k')
scatter(s(end,2)/1e3,T(end,2),msz,'or','MarkerFaceColor','r','MarkerEdgeColor','k')

px = [s(:,1)/1e3; s(end,2)/1e3];
py = [T(:,1); T(end,2)];
for i=1:length(px)
    str2 = sprintf('%i',i);
    if i==2
        px(i) = px(i)*30;
        py(i) = py(i)*1.3;
    end
    text(px(i)-.4, py(i), str2, 'Interpreter', 'latex')
end

pbaspect([1.2 1 1])
set(gca,'LineWidth',alw)
xlim([-1 6])
xlabel('$s~(\mathrm{kJ/kg\cdot K})$')
ylim([0 700])
ylab=ylabel('$T~(\mathrm{K})$');
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

% print('Ts_baseline','-djpeg','-r400')


% fs=15;
% n = length(p);
% h=text(s(1,1)/1e3,T(1,1),'$20~\mathrm{MPa}$','FontSize',fs);
% set(h,'Rotation',75);
% h=text(3,400,'$3~\mathrm{MPa}$','FontSize',fs);
% set(h,'Rotation',75);
% h=text(3.3,270,'$1.5~\mathrm{MPa}$','FontSize',fs);
% set(h,'Rotation',75);
% h=text(3.75,380,'$0.2~\mathrm{MPa}$','FontSize',fs);
% set(h,'Rotation',75);
% h=text(4.35,380,'$0.1~\mathrm{MPa}$','FontSize',fs);
% set(h,'Rotation',75);
% h=text(-0.5,450,'AA');
% px = data(1:n,1)/1e3;
% py = data(1:n,2);
% for i=1:length(px)
%     str2 = sprintf('%i',i);
%     text(px(i)-.4, py(i), str2, 'Interpreter', 'latex')
% end
