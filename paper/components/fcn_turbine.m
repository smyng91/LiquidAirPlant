function [W, T_out, X_dest] = fcn_turbine( TURBINE, AMBIENT )

% turbine model

try 
    
h_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',TURBINE.p_in,'T',TURBINE.T_in,TURBINE.FLUID);
s_in = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE.p_in,'T',TURBINE.T_in,TURBINE.FLUID);
h_out_is = py.CoolProp.CoolProp.PropsSI('Hmass','P',TURBINE.p_out,'Smass',s_in,TURBINE.FLUID);
h_out = h_in-TURBINE.eta_is*(h_in-h_out_is);
s_out = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE.p_out,'Hmass',h_out,TURBINE.FLUID);

% output
W = TURBINE.eta_is*TURBINE.mdot*(h_in-h_out_is);
T_out = py.CoolProp.CoolProp.PropsSI('T','P',TURBINE.p_out,'Hmass',h_out,TURBINE.FLUID);
X_dest = TURBINE.mdot*AMBIENT.T0*(s_out-s_in);

catch
    
    T_out = 0;
    X_dest = 0;
    W = 0;
    
end

end
