function [y, PLANT] = model_AA( T, PARAM )
 global data_Ts_AA
%
% Liquid air plant model
% Configuration - Ambient Air
%
options = optimoptions('fsolve','Display','none');

PREHEATER.T_c_out = T(1)*PARAM.T0;
PREHEATER.T_h_out = T(2)*PARAM.T0;
REGENERATOR.T_c_out = T(3)*PARAM.T0;
REGENERATOR.T_h_out = T(4)*PARAM.T0;

% pump
PUMP.FLUID = PARAM.FLUID;
PUMP.p_in = PARAM.p_store;
PUMP.T_in = PARAM.T_store;
PUMP.mdot = PARAM.mdot_la;
PUMP.eta_is = PARAM.eta_pump;
PUMP.p_out = PARAM.p_pump;
[PUMP.W, PUMP.T_out, PUMP.X_dest] = fcn_pump( PUMP, PARAM );

% preheater
PREHEATER.FLUID_h = PARAM.FLUID;
PREHEATER.FLUID_c = PARAM.FLUID;
PREHEATER.mdot_c = PUMP.mdot;
PREHEATER.mdot_h = PARAM.mdot_preheater;
PREHEATER.T_h_in = PARAM.T0;
PREHEATER.T_c_in = PUMP.T_out;
PREHEATER.p_c_in = PUMP.p_out;
PREHEATER.p_h_in = PARAM.p0;
PREHEATER.p_c_out = PREHEATER.p_c_in;
PREHEATER.p_h_out = PREHEATER.p_h_in;
PREHEATER.NTU = PARAM.NTU_pre;
[x, PREHEATER.X_dest] = fcn_hx( PREHEATER, PARAM );
y(1) = abs( x(1) - PREHEATER.T_c_out )/PARAM.T0;
y(2) = abs( x(2) - PREHEATER.T_h_out )/PARAM.T0;

TURBINE_HP.FLUID = PARAM.FLUID;
TURBINE_HP.eta_is = PARAM.eta_turbine;
TURBINE_HP.p_in = PREHEATER.p_c_out;
TURBINE_HP.p_out = 3e6;
TURBINE_HP.T_in = PREHEATER.T_c_out;
TURBINE_HP.mdot = PREHEATER.mdot_c;
[TURBINE_HP.W, TURBINE_HP.T_out, TURBINE_HP.X_dest] = fcn_turbine( TURBINE_HP, PARAM );

REGENERATOR.FLUID_h = PARAM.FLUID;
REGENERATOR.FLUID_c = PARAM.FLUID;
REGENERATOR.mdot_c = TURBINE_HP.mdot;
REGENERATOR.mdot_h = PARAM.mdot_regen;
REGENERATOR.T_h_in = PARAM.T0;
REGENERATOR.T_c_in = TURBINE_HP.T_out;
REGENERATOR.p_c_in = TURBINE_HP.p_out;
REGENERATOR.p_h_in = PARAM.p0;
REGENERATOR.p_c_out = REGENERATOR.p_c_in;
REGENERATOR.p_h_out = REGENERATOR.p_h_in;
REGENERATOR.NTU = PARAM.NTU_reg;
[x, REGENERATOR.X_dest] = fcn_hx( REGENERATOR, PARAM );
y(3) = abs( x(1) - REGENERATOR.T_c_out )/PARAM.T0;
y(4) = abs( x(2) - REGENERATOR.T_h_out )/PARAM.T0;

TURBINE_LP.FLUID = PARAM.FLUID;
TURBINE_LP.eta_is = PARAM.eta_turbine;
TURBINE_LP.p_in = REGENERATOR.p_c_out;
TURBINE_LP.p_out = PARAM.p0;
TURBINE_LP.T_in = REGENERATOR.T_c_out;
TURBINE_LP.mdot = REGENERATOR.mdot_c;
[TURBINE_LP.W, TURBINE_LP.T_out, TURBINE_LP.X_dest] = fcn_turbine( TURBINE_LP, PARAM );

PLANT.W_t = TURBINE_HP.W + TURBINE_LP.W;
PLANT.W_net = TURBINE_HP.W + TURBINE_LP.W - PUMP.W;
PLANT.W_in = (0.5*1.e3*3600)*PUMP.mdot;
PLANT.eta_rt = PLANT.W_net/PLANT.W_in;

try 
    s1 = py.CoolProp.CoolProp.PropsSI('Smass','P',PUMP.p_in,'T',PUMP.T_in,PUMP.FLUID);
    s2 = py.CoolProp.CoolProp.PropsSI('Smass','P',PUMP.p_out,'T',PUMP.T_out,PUMP.FLUID);
    s3 = py.CoolProp.CoolProp.PropsSI('Smass','P',PREHEATER.p_c_out,'T',PREHEATER.T_c_out,PREHEATER.FLUID_c);
    s4 = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE_HP.p_out,'T',TURBINE_HP.T_out,TURBINE_HP.FLUID);
    s5 = py.CoolProp.CoolProp.PropsSI('Smass','P',REGENERATOR.p_c_out,'T',REGENERATOR.T_c_out,REGENERATOR.FLUID_c);
    s6 = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE_LP.p_out,'T',TURBINE_LP.T_out,TURBINE_LP.FLUID);
    data_Ts_AA = [s1,PUMP.T_in; s2,PUMP.T_out; s3,PREHEATER.T_c_out; s4,TURBINE_HP.T_out; s5,REGENERATOR.T_c_out; s6,TURBINE_LP.T_out];
catch
    y = 1000*ones(length(T),1);
end
    