function [W, T_out, X_dest] = fcn_pump( PUMP, AMBIENT )

% pump model

% error/exception handling is needed to avoid divergence in physically
% infeasible design spaces, where CoolProp will output error and stop the
% code.
try
    
h_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',PUMP.p_in,'T',PUMP.T_in,PUMP.FLUID);
s_in = py.CoolProp.CoolProp.PropsSI('Smass','P',PUMP.p_in,'T',PUMP.T_in,PUMP.FLUID);
h_out_is = py.CoolProp.CoolProp.PropsSI('Hmass','P',PUMP.p_out,'Smass',s_in,PUMP.FLUID);
h_out = h_in+1/PUMP.eta_is*(h_out_is-h_in);
s_out = py.CoolProp.CoolProp.PropsSI('Smass','P',PUMP.p_out,'Hmass',h_out,PUMP.FLUID);

% output
W = PUMP.mdot/PUMP.eta_is * (h_out_is - h_in);
T_out = py.CoolProp.CoolProp.PropsSI('T','P',PUMP.p_out,'Hmass',h_out,PUMP.FLUID);
X_dest = PUMP.mdot*AMBIENT.T0*(s_out-s_in);

catch
    
    W = 0;
    T_out = 0;
    X_dest = 0;
    
end

