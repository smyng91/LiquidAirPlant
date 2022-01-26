function [dT, p_out, X_sun, X_dest, eta_c] = fcn_PTC_validation( T, PTC, AMBIENT )

% PTC model

PTC.T_r = T(1);
PTC.T_c = T(2);
PTC.T_out = T(3);

% try
    
    % evaluate properties
    rho_air = py.CoolProp.CoolProp.PropsSI('Dmass','P',AMBIENT.p0,'T',AMBIENT.T0,AMBIENT.FLUID);
    mu_air = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',AMBIENT.p0,'T',AMBIENT.T0,AMBIENT.FLUID);
    k_air = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','P',AMBIENT.p0,'T',AMBIENT.T0,AMBIENT.FLUID);
    Re_air = rho_air*AMBIENT.u0*PTC.d_c/mu_air;
    A_r_cs = pi*(PTC.d_r/1.05)^2/4;
    
    % convection, ambient
    if Re_air < 1e3
        Nu = 0.4+0.54*Re_air^0.52;
    else
        Nu = 0.3*Re_air^0.6;
    end
    h_w = Nu*k_air/PTC.d_c;
    
    % radiation, ambient
    h_r_c0 = PTC.e_c*PTC.sigma*(PTC.T_c+AMBIENT.T0)*(PTC.T_c^2+AMBIENT.T0^2);
    
    % radiation, absorber-cover
    h_r_rc = PTC.sigma*(PTC.T_r^2+PTC.T_c^2)*(PTC.T_r+PTC.T_c)/(1/PTC.e_r+PTC.A_r/PTC.A_c*(1/PTC.e_c-1));
    
    T_mean = (PTC.T_in+PTC.T_out)/2;
    rho_HTF = py.CoolProp.CoolProp.PropsSI('Dmass','P',PTC.p_in,'T',T_mean,PTC.FLUID);
    mu_HTF = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',PTC.p_in,'T',T_mean,PTC.FLUID);
    cp_HTF = py.CoolProp.CoolProp.PropsSI('CPMASS','P',PTC.p_in,'T',T_mean,PTC.FLUID);
    k_HTF = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','P',PTC.p_in,'T',T_mean,PTC.FLUID);
    Pr_HTF = py.CoolProp.CoolProp.PropsSI('PRANDTL','P',PTC.p_in,'T',T_mean,PTC.FLUID);
    u_HTF = PTC.mdot/(rho_HTF*A_r_cs);
    Re_HTF = rho_HTF*u_HTF*PTC.d_r/1.05/mu_HTF;
    
    % convection, HTF-absorber
    if Re_HTF > 2300
        Nu = 0.023*(Re_HTF)^0.8*Pr_HTF^0.4;
    else
        Nu = 4.364;
    end
    h_fi = Nu*k_HTF/PTC.d_r/1.05;
    
    % overall heat loss
    k = 280;
    U_l = 1/( PTC.A_r/((h_w+h_r_c0)*PTC.A_c) + 1/h_r_rc );
    U_o = 1/( 1/U_l+1.05/h_fi+PTC.d_r*log(1.05)/(2*k) );
    
    Fprime = U_o/U_l;
    Fr = PTC.mdot*cp_HTF/(PTC.A_r*U_l)*(1-exp(-U_l*Fprime*PTC.A_r/(PTC.mdot*cp_HTF)));
    Qu = Fr*(AMBIENT.Gb*PTC.eta_opt*PTC.A_ap-PTC.A_r*U_l*(PTC.T_in-AMBIENT.T0));
    dT(1) = abs(PTC.T_in + Qu/(PTC.mdot*cp_HTF)-PTC.T_out);
    dT(2) = abs( (PTC.A_r*h_r_rc*PTC.T_r+PTC.A_c*(h_r_c0+h_w)*AMBIENT.T0)/(PTC.A_r*h_r_rc+PTC.A_c*(h_r_c0+h_w))-PTC.T_c);
    dT(3) = abs((AMBIENT.Gb*PTC.eta_opt*PTC.A_ap-Qu)/(PTC.A_r*U_l)+AMBIENT.T0 - PTC.T_r);
    
    % pumping power
    % Haaland equation for friction
    rough = 1e-3;
    f = ((-1.8*log((rough/PTC.d_r/1.05/3.7)^1.11+6.9/Re_HTF))^-1)^2;
    dp = PTC.L*f*rho_HTF/2*u_HTF^2/PTC.d_r/1.05;
    p_out = PTC.p_in - dp;
    
    p_mean = (PTC.p_in+p_out)/2;
    
    % exergy
    s_in = py.CoolProp.CoolProp.PropsSI('Smass','P',PTC.p_in,'T',PTC.T_in,PTC.FLUID);
    s_out = py.CoolProp.CoolProp.PropsSI('Smass','P',p_out,'T',PTC.T_out,PTC.FLUID);
    rho_f = py.CoolProp.CoolProp.PropsSI('Dmass','P',p_mean,'T',T_mean,PTC.FLUID);
    X_sun = AMBIENT.Gb*PTC.A_ap*( 1-4/3*(AMBIENT.T0/AMBIENT.T_sun)+1/3*(AMBIENT.T0/AMBIENT.T_sun)^4 );
    X_dest_dP = PTC.mdot*AMBIENT.T0*dp/rho_f*log(PTC.T_out/PTC.T_in)/(PTC.T_out-PTC.T_in);
    X_dest_Q1 = X_sun*PTC.eta_opt-AMBIENT.Gb*PTC.A_ap*PTC.eta_opt*(1-AMBIENT.T0/PTC.T_r);
    X_dest_Q2 = PTC.mdot*AMBIENT.T0*(s_out-s_in)-AMBIENT.T0*Qu/PTC.T_r;
    X_dest = X_dest_Q1 + X_dest_Q2 + X_dest_dP;
    eta_c = Fr*(PTC.eta_opt - U_l*(PTC.T_in-AMBIENT.T0)/AMBIENT.Gb/PTC.Cr);
%     
% catch
%     
%     dT = 1000*ones(3,1);
%     X_dest = 0;
%     X_sun = 0;
%     p_out = PTC.p_in;
%     eta_c = 0;
%     
% end
% 
