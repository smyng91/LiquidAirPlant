clc; clear; close all; clear fcn;
alw = 1.2;    % AxesLineWidth
fsz = 20;      % Fontsize
ftype = 'Times New Roman';  % Font type
lw = 1.5;      % LineWidth
msz = 15;       % MarkerSize
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

% PTC
PTC.FLUID = 'INCOMP::S800';
PTC.eta_opt = 0.726;
PTC.p_in = 0.2e6;
PTC.d_c = 0.105;
PTC.d_r = 0.066;
PTC.e_c = 0.9;
PTC.e_r = 0.14;
PTC.sigma = 5.67e-8;
PTC.rim = 50;
PTC.L = 8;
PTC.L_f = 1.84;
PTC.A_ap = 4*PTC.L_f*tand(PTC.rim/2)*PTC.L;
PTC.A_r = pi()*PTC.d_r*PTC.L;
PTC.A_c = pi()*PTC.d_c*PTC.L;
PTC.Cr = PTC.A_ap/PTC.A_r;
AMBIENT.p0 = 101325*0.83;
AMBIENT.FLUID = 'AIR';
AMBIENT.T_sun = 5770;

%define vectors to match experimental parameters for efficiency in vacuum
DNI_vec=[933.7;968.2;982.3;909.5;937.9;880.6;903.2;920.9]; %direct normal insolation (W/m^2)
U_ext_vec = [2.6;3.7;2.5;3.3;1.0;2.9;4.2;2.6;]; %wind speed (m/s)
T_ext_vec = [21.2;22.4;24.3;26.2;28.8;27.5;31.1;29.5]+273.15; %ambient temperature (m/s)
T_fi_vec = [102.2;151;197.5;250.7;297.8;299;355.9;379.5]+273.15; %inlet temperature (K)
DeltaT_eff_vac = [91.9;139.8;184.3;233.9;278.6;280.7;334.1;359.4]; %temp above ambient (K)
Vdot_f_vec = [47.7;47.8;49.1;54.7;55.5;55.6;56.3;56.8]*10; %volumetric flow rate (L/min)
Eff_Vac_Dudley = [72.51;70.9;70.17;70.25;67.98;68.92;63.82;62.34]; %measured efficiency
Eff_Vac_Dudley_Error = [1.95;1.92;1.81;1.90;1.86;2.06;2.36;2.41]; %associated error

x0 = [400,300,340];
for i=1:length(DNI_vec)
    PTC.T_in = T_fi_vec(i);
    AMBIENT.Gb = DNI_vec(i);
    AMBIENT.u0 = U_ext_vec(i);
    AMBIENT.T0 = T_ext_vec(i);
    AMBIENT.h0 = py.CoolProp.CoolProp.PropsSI('Hmass','P',AMBIENT.p0,'T',AMBIENT.T0, AMBIENT.FLUID);
    AMBIENT.s0 = py.CoolProp.CoolProp.PropsSI('Smass','P',AMBIENT.p0,'T',AMBIENT.T0, AMBIENT.FLUID);
    rho = py.CoolProp.CoolProp.PropsSI('Dmass','P',PTC.p_in,'T',PTC.T_in,'INCOMP::S800');
    PTC.mdot = rho*Vdot_f_vec(i)*1.667*10^-5;
    func = @(x) fcn_PTC_validation( x, PTC, AMBIENT );
    [T, fval, exitflag, output] = fsolve(func, x0);

    PTC.T_r = T(1);
    PTC.T_c = T(2);
    PTC.T_out = T(3);
    [T, PTC.p_out, PTC.X_sun, PTC.X_dest, PTC.eta_c] = fcn_PTC( PTC, AMBIENT );
    eta_c(i) = PTC.eta_c;
end




figure('units','normalized','position',[0.2 0.2 .5 .5])

h2 = errorbar(DeltaT_eff_vac,Eff_Vac_Dudley,Eff_Vac_Dudley_Error,'o','Color','k',...
    'LineWidth',2,'MarkerSize',msz,'MarkerFaceColor','w')
hold on
ylim([55 80])
grid minor
h1 = plot(DeltaT_eff_vac,eta_c*100,'o','Color','k','MarkerSize',msz,'LineWidth',2,'MarkerFaceColor','r')
ylab = ylabel('$\mathrm{\eta_{c}}$','Interpreter','LaTex');
pbaspect([1.2 1 1])
set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5]);
xlabel('$\Delta T~\mathrm{(K)}$')
legend([h1,h2],'Present model','Dudley et al.')

% get handle to current axes
a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'XColor','k');
% linkaxes([a b])
pbaspect([1.2 1 1])
axes(a)

