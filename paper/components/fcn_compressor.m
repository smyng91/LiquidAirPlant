function [W, T_out, X_dest] = fcn_compressor( COMPRESSOR, AMBIENT )

% compressor model

h_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',COMPRESSOR.p_in,'T',COMPRESSOR.T_in,COMPRESSOR.FLUID);
s_in = py.CoolProp.CoolProp.PropsSI('Smass','P',COMPRESSOR.p_in,'T',COMPRESSOR.T_in,COMPRESSOR.FLUID);
h_out_is = py.CoolProp.CoolProp.PropsSI('Hmass','P',COMPRESSOR.p_out,'Smass',s_in,COMPRESSOR.FLUID);
h_out = h_in+1/COMPRESSOR.eta_is*(h_out_is-h_in);
s_out = py.CoolProp.CoolProp.PropsSI('Smass','P',COMPRESSOR.p_out,'Hmass',h_out,COMPRESSOR.FLUID);

% output
W = COMPRESSOR.mdot/COMPRESSOR.eta_is * (h_out_is - h_in);
T_out = py.CoolProp.CoolProp.PropsSI('T','P',COMPRESSOR.p_out,'Hmass',h_out,COMPRESSOR.FLUID);
X_dest = COMPRESSOR.mdot*AMBIENT.T0*(s_out-s_in);

end