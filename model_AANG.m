function [y, PLANT ] = model_AANG( T, PARAM )
% LAPP model - Ambient Air/Natural Gas

global data_Ts_AANG

PREHEATER.T_c_out = T(1)*PARAM.T0;
PREHEATER.T_h_out = T(2)*PARAM.T0;
REGENERATOR.T_c_out = T(3)*PARAM.T0;
REGENERATOR.T_h_out = T(4)*PARAM.T0;
TURBINE_LP.T_out = T(5)*PARAM.T0;

% error/exception handling is needed to avoid divergence in physically
% infeasible design spaces, where CoolProp will output error and stop the
% code.
try
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

    % regenerator
    REGENERATOR.FLUID_h = PARAM.FLUID;
    REGENERATOR.FLUID_c = PARAM.FLUID;
    REGENERATOR.mdot_c = PREHEATER.mdot_c;
    REGENERATOR.mdot_h = REGENERATOR.mdot_c;
    REGENERATOR.T_h_in = TURBINE_LP.T_out;
    REGENERATOR.T_c_in = PREHEATER.T_c_out;
    REGENERATOR.p_c_in = PREHEATER.p_c_out;
    REGENERATOR.p_h_in = PARAM.p0;
    REGENERATOR.p_c_out = REGENERATOR.p_c_in;
    REGENERATOR.p_h_out = REGENERATOR.p_h_in;
    REGENERATOR.NTU = PARAM.NTU_reg;
    [x, REGENERATOR.X_dest] = fcn_hx( REGENERATOR, PARAM );
    y(3) = abs( x(1) - REGENERATOR.T_c_out )/PARAM.T0;
    y(4) = abs( x(2) - REGENERATOR.T_h_out )/PARAM.T0;

    % Combustor 1
    COMBUST1.FLUID = PARAM.FLUID;
    COMBUST1.T_in = REGENERATOR.T_c_out;
    COMBUST1.p_in = REGENERATOR.p_c_out;
    COMBUST1.p_out = COMBUST1.p_in;
    COMBUST1.mdot = REGENERATOR.mdot_c;
    COMBUST1.Q = PARAM.Q_NG/2;
    COMBUST1.NTU = PARAM.NTU_combus1;
    [ COMBUST1.T_out, COMBUST1.X_dest ] = fcn_combustor( COMBUST1, PARAM );

    % HP turbine
    TURBINE_HP.FLUID = PARAM.FLUID;
    TURBINE_HP.eta_is = PARAM.eta_turbine;
    TURBINE_HP.p_in = COMBUST1.p_out;
    TURBINE_HP.p_out = 3e6;
    TURBINE_HP.T_in = COMBUST1.T_out;
    TURBINE_HP.mdot = COMBUST1.mdot;
    [TURBINE_HP.W, TURBINE_HP.T_out, TURBINE_HP.X_dest] = fcn_turbine( TURBINE_HP, PARAM );

    % Combustor 2
    COMBUST2.FLUID = PARAM.FLUID;
    COMBUST2.T_in = TURBINE_HP.T_out;
    COMBUST2.p_in = TURBINE_HP.p_out;
    COMBUST2.p_out = COMBUST2.p_in;
    COMBUST2.mdot = TURBINE_HP.mdot;
    COMBUST2.Q = PARAM.Q_NG/2;
    COMBUST2.NTU = PARAM.NTU_combus2;
    [ COMBUST2.T_out, COMBUST2.X_dest ] = fcn_combustor( COMBUST2, PARAM );

    % LP turbine
    TURBINE_LP.FLUID = PARAM.FLUID;
    TURBINE_LP.eta_is = PARAM.eta_turbine;
    TURBINE_LP.p_in = COMBUST2.p_out;
    TURBINE_LP.p_out = PARAM.p0;
    TURBINE_LP.T_in = COMBUST2.T_out;
    TURBINE_LP.mdot = COMBUST2.mdot;
    [TURBINE_LP.W, TURBINE_LP.T_out, TURBINE_LP.X_dest] = fcn_turbine( TURBINE_LP, PARAM );
    y(5) = abs( REGENERATOR.T_h_in - TURBINE_LP.T_out )/PARAM.T0;

    % Equations from M. Antonelli et al. / Applied Energy 194 (2017)
    % 522–529,
    PLANT.W_t = TURBINE_HP.W + TURBINE_LP.W;
    PLANT.W_net = TURBINE_HP.W + TURBINE_LP.W - PUMP.W;
    PLANT.W_in = PARAM.e_liq*PUMP.mdot + PARAM.eta_cc*PARAM.Q_NG;
    PLANT.eta_rt = PLANT.W_net/PLANT.W_in;

    s1 = py.CoolProp.CoolProp.PropsSI('Smass','P',PUMP.p_in,'T',PUMP.T_in,PUMP.FLUID);
    s2 = py.CoolProp.CoolProp.PropsSI('Smass','P',PUMP.p_out,'T',PUMP.T_out,PUMP.FLUID);
    s3 = py.CoolProp.CoolProp.PropsSI('Smass','P',PREHEATER.p_c_out,'T',PREHEATER.T_c_out,PREHEATER.FLUID_c);
    s4 = py.CoolProp.CoolProp.PropsSI('Smass','P',REGENERATOR.p_c_out,'T',REGENERATOR.T_c_out,REGENERATOR.FLUID_c);
    s5 = py.CoolProp.CoolProp.PropsSI('Smass','P',COMBUST1.p_out,'T',COMBUST1.T_out,COMBUST1.FLUID);
    s6 = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE_HP.p_out,'T',TURBINE_HP.T_out,TURBINE_HP.FLUID);
    s7 = py.CoolProp.CoolProp.PropsSI('Smass','P',COMBUST2.p_out,'T',COMBUST2.T_out,COMBUST2.FLUID);
    s8 = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE_LP.p_out,'T',TURBINE_LP.T_out,TURBINE_LP.FLUID);
    s9 = py.CoolProp.CoolProp.PropsSI('Smass','P',REGENERATOR.p_h_out,'T',REGENERATOR.T_h_out,REGENERATOR.FLUID_h);
    data_Ts_AANG = [...
        s1,PUMP.T_in;
        s2,PUMP.T_out;
        s3,PREHEATER.T_c_out;
        s4,REGENERATOR.T_c_out;
        s5,COMBUST1.T_out;
        s6,TURBINE_HP.T_out;
        s7,COMBUST2.T_out;
        s8,TURBINE_LP.T_out;
        s9,REGENERATOR.T_h_out];

catch

    y = ones(length(T),1);

end

end


