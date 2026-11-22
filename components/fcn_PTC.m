function [T, p_out, X_sun, X_dest, eta_c, Qu, dp] = fcn_PTC( PTC, AMBIENT )
%FCN_PTC Parabolic trough collector (single receiver tube with glass cover).
%
%   [T, p_out, X_sun, X_dest, eta_c, Qu, dp] = FCN_PTC( PTC, AMBIENT )
%
%   PTC fields: FLUID, mdot [kg/s], p_in [Pa], T_in [K], eta_opt [-],
%               d_c, d_r [m] (cover / receiver outer diameters),
%               d_ri [m] (receiver inner diameter, default d_r/1.05),
%               k_r [W/m/K] (receiver wall conductivity, default 280),
%               roughness [m] (receiver inner wall, default 1e-3),
%               e_c, e_r [-], sigma, L [m], A_ap, A_r, A_c [m2], Cr [-],
%               and the solver's current guesses T_r, T_c, T_out [K]
%   T       computed [T_r, T_c, T_out] [K]
%   p_out   outlet pressure [Pa]
%   X_sun   exergy of the incident solar radiation [W]
%   X_dest  exergy destruction [W]
%   eta_c   collector thermal efficiency [-]
%   Qu      useful heat gain [W]
%   dp      receiver pressure drop [Pa]

d_r = PTC.d_r;
d_ri = field_or( PTC, 'd_ri', d_r/1.05 );
k_r = field_or( PTC, 'k_r', 280 );
roughness = field_or( PTC, 'roughness', 1e-3 );
T0 = AMBIENT.T0;
S = AMBIENT.Gb*PTC.eta_opt*PTC.A_ap;   % absorbed solar radiation [W]

% Petela exergy of solar radiation
X_sun = AMBIENT.Gb*PTC.A_ap*( 1 - 4/3*(T0/AMBIENT.T_sun) + 1/3*(T0/AMBIENT.T_sun)^4 );

% CoolProp failures (physically infeasible states) are caught so that the
% solver can recover; any other error is a genuine bug and is rethrown.
try
    % convection, cover-ambient (wind across the glass cover)
    [rho_air, mu_air, k_air] = fcn_props( AMBIENT.FLUID, 'PT', AMBIENT.p0, T0, 'D', 'MU', 'K' );
    Re_air = rho_air*AMBIENT.u0*PTC.d_c/mu_air;
    if Re_air < 1e3
        Nu = 0.4 + 0.54*Re_air^0.52;
    else
        Nu = 0.3*Re_air^0.6;
    end
    h_w = Nu*k_air/PTC.d_c;

    % radiation, cover-sky
    h_r_c0 = PTC.e_c*PTC.sigma*(PTC.T_c + T0)*(PTC.T_c^2 + T0^2);

    % radiation, receiver-cover
    h_r_rc = PTC.sigma*(PTC.T_r^2 + PTC.T_c^2)*(PTC.T_r + PTC.T_c) / ...
        ( 1/PTC.e_r + PTC.A_r/PTC.A_c*(1/PTC.e_c - 1) );

    % heat transfer fluid, evaluated at the mean temperature
    T_mean = (PTC.T_in + PTC.T_out)/2;
    [rho_HTF, mu_HTF, cp_HTF, k_HTF, Pr_HTF] = ...
        fcn_props( PTC.FLUID, 'PT', PTC.p_in, T_mean, 'D', 'MU', 'CP', 'K', 'PR' );
    u_HTF = PTC.mdot/(rho_HTF*pi*d_ri^2/4);
    Re_HTF = rho_HTF*u_HTF*d_ri/mu_HTF;

    % convection, HTF-receiver (Dittus-Boelter / fully developed laminar)
    if Re_HTF > 2300
        Nu = 0.023*Re_HTF^0.8*Pr_HTF^0.4;
    else
        Nu = 4.364;
    end
    h_fi = Nu*k_HTF/d_ri;

    % overall heat loss and heat transfer coefficients (receiver outer area)
    U_l = 1/( PTC.A_r/((h_w + h_r_c0)*PTC.A_c) + 1/h_r_rc );
    U_o = 1/( 1/U_l + d_r/(d_ri*h_fi) + d_r*log(d_r/d_ri)/(2*k_r) );

    % collector efficiency and heat removal factors, useful heat gain
    Fprime = U_o/U_l;
    mcp = PTC.mdot*cp_HTF;
    Fr = mcp/(PTC.A_r*U_l)*( 1 - exp(-U_l*Fprime*PTC.A_r/mcp) );
    Qu = Fr*( S - PTC.A_r*U_l*(PTC.T_in - T0) );
    T_out = PTC.T_in + Qu/mcp;
    T_c = (PTC.A_r*h_r_rc*PTC.T_r + PTC.A_c*(h_r_c0 + h_w)*T0) / ...
        (PTC.A_r*h_r_rc + PTC.A_c*(h_r_c0 + h_w));
    T_r = (S - Qu)/(PTC.A_r*U_l) + T0;
    T = [T_r, T_c, T_out];
    eta_c = Fr*( PTC.eta_opt - U_l*(PTC.T_in - T0)/AMBIENT.Gb/PTC.Cr );

    % pressure drop: Darcy friction factor from Haaland (turbulent) or
    % 64/Re (laminar), based on the receiver inner diameter
    if Re_HTF > 2300
        f = ( -1.8*log10( (roughness/d_ri/3.7)^1.11 + 6.9/Re_HTF ) )^-2;
    else
        f = 64/Re_HTF;
    end
    dp = f*PTC.L/d_ri*rho_HTF*u_HTF^2/2;
    p_out = PTC.p_in - dp;

    % exergy destruction: sun -> receiver, receiver -> HTF, friction
    s_in = fcn_props( PTC.FLUID, 'PT', PTC.p_in, PTC.T_in, 'S' );
    s_out = fcn_props( PTC.FLUID, 'PT', p_out, T_out, 'S' );
    if abs(T_out - PTC.T_in) > 1e-6*PTC.T_in
        inv_T_lm = log(T_out/PTC.T_in)/(T_out - PTC.T_in);
    else
        inv_T_lm = 1/PTC.T_in;
    end
    X_dest_dP = PTC.mdot*T0*dp/rho_HTF*inv_T_lm;
    X_dest_Q1 = X_sun*PTC.eta_opt - S*(1 - T0/T_r);
    X_dest_Q2 = PTC.mdot*T0*(s_out - s_in) - T0*Qu/T_r;
    X_dest = X_dest_Q1 + X_dest_Q2 + X_dest_dP;
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    T = [1000, 1000, 1000];
    p_out = PTC.p_in;
    X_dest = NaN;
    eta_c = NaN;
    Qu = NaN;
    dp = NaN;
end

end

function v = field_or( S, name, default )
if isfield(S, name)
    v = S.(name);
else
    v = default;
end
end
