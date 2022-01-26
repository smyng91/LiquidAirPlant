function [T_out, X_dest] = fcn_combustor( COMBUST, AMBIENT )

% combustor model

% error/exception handling is needed to avoid divergence in physically
% infeasible design spaces, where CoolProp will output error and stop the
% code.
try
    cp = py.CoolProp.CoolProp.PropsSI('CPMASS','P',COMBUST.p_in,'T',COMBUST.T_in,COMBUST.FLUID);
    e = 1-exp(-COMBUST.NTU);
    T_out = COMBUST.T_in + e*COMBUST.Q/(COMBUST.mdot*cp);
    Tavg = (COMBUST.T_in + T_out)/2;

    h_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',COMBUST.p_in,'T',COMBUST.T_in,COMBUST.FLUID);
    h_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',COMBUST.p_out,'T',T_out,COMBUST.FLUID);
    s_in = py.CoolProp.CoolProp.PropsSI('Smass','P',COMBUST.p_in,'T',COMBUST.T_in,COMBUST.FLUID);
    s_out = py.CoolProp.CoolProp.PropsSI('Smass','P',COMBUST.p_out,'T',T_out,COMBUST.FLUID);
    X_dest = AMBIENT.T0*COMBUST.mdot*( h_in-h_out - AMBIENT.T0*(s_in-s_out) ) + COMBUST.Q*(1-AMBIENT.T0/Tavg);
    
catch
    X_dest = 0;
    T_out = 0;
end

end



