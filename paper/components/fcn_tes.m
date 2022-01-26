function [dT, X_dest] = fcn_tes( TES, AMBIENT )

% TES model

try
    
    h_ch_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',TES.p_ch_in,'T',TES.T_ch_in,TES.FLUID_ch);
    h_ch_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',TES.p_ch_out,'T',TES.T,TES.FLUID_ch);
    h_dch_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',TES.p_dch_in,'T',TES.T_dch_in,TES.FLUID_dch);
    h_dch_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',TES.p_dch_out,'T',TES.T,TES.FLUID_dch);
    c =  py.CoolProp.CoolProp.PropsSI('CPMASS ','P',TES.p,'T',TES.T,TES.FLUID_ch);
    
    s_ch_in = py.CoolProp.CoolProp.PropsSI('Smass','P',TES.p_ch_in,'T',TES.T_ch_in,TES.FLUID_ch);
    s_ch_out = py.CoolProp.CoolProp.PropsSI('Smass','P',TES.p_ch_out,'T',TES.T,TES.FLUID_ch);
    s_dch_in = py.CoolProp.CoolProp.PropsSI('Smass','P',TES.p_dch_in,'T',TES.T_dch_in,TES.FLUID_dch);
    s_dch_out = py.CoolProp.CoolProp.PropsSI('Smass','P',TES.p_dch_out,'T',TES.T,TES.FLUID_dch);
    
    dT = 1/TES.M*(TES.mdot_ch*( h_ch_in-h_ch_out)+TES.mdot_dch*(h_dch_in-h_dch_out) );
    X_dest = AMBIENT.T0*( (TES.M*c/TES.T)*dT + TES.mdot_ch*( s_ch_out-s_ch_in ) + TES.mdot_dch*(s_dch_in-s_dch_out)   );
    
catch
    
    dT = 0;
    X_dest = 0;
    
end

end