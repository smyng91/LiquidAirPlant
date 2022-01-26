function [y, PLANT ] = model_RAANG( T, PARAM )
% LAPP model - Recuperative Ambient Air/Natural Gas

global data_Ts_RAANG

PREHEATER.T_c_out = T(1)*PARAM.T0;
PREHEATER.T_h_out = T(2)*PARAM.T0;
REGENERATOR.T_c_out = T(3)*PARAM.T0;
REGENERATOR.T_h_out = T(4)*PARAM.T0;
RECUPERATOR.T_c_out = T(5)*PARAM.T0;
RECUPERATOR.T_h_out = T(6)*PARAM.T0;
TURBINE_LP.T_out = T(7)*PARAM.T0;
TURBINE_RECUP.T_out = T(8)*PARAM.T0;

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
    COMBUST1.Q = PARAM.Q_NG(1)*2/3;
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
    COMBUST2.Q = PARAM.Q_NG(2)/3;
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
    y(7) =  abs( REGENERATOR.T_h_in - TURBINE_LP.T_out)/PARAM.T0;
    
    % RECUPERATIVE CYCLE
    
    % compressor
    COMPRESSOR.FLUID = PARAM.FLUID;
    COMPRESSOR.T_in = PREHEATER.T_h_out;
    COMPRESSOR.p_in = PREHEATER.p_h_out;
    COMPRESSOR.p_out = 1.5e6;
    COMPRESSOR.eta_is = PARAM.eta_comp;
    COMPRESSOR.mdot = PREHEATER.mdot_h;
    [COMPRESSOR.W, COMPRESSOR.T_out, COMPRESSOR.X_dest] = fcn_compressor( COMPRESSOR, PARAM );
    
    % recuperator
    RECUPERATOR.FLUID_h = PARAM.FLUID;
    RECUPERATOR.FLUID_c = PARAM.FLUID;
    RECUPERATOR.mdot_c = COMPRESSOR.mdot;
    RECUPERATOR.mdot_h = RECUPERATOR.mdot_c;
    RECUPERATOR.T_h_in = TURBINE_RECUP.T_out;
    RECUPERATOR.T_c_in = COMPRESSOR.T_out;
    RECUPERATOR.p_c_in = COMPRESSOR.p_out;
    RECUPERATOR.p_h_in = PARAM.p0;
    RECUPERATOR.p_c_out = RECUPERATOR.p_c_in;
    RECUPERATOR.p_h_out = RECUPERATOR.p_h_in;
    RECUPERATOR.NTU = PARAM.NTU_rec;
    [x, RECUPERATOR.X_dest] = fcn_hx( RECUPERATOR, PARAM );
    y(5) = abs( x(1) - RECUPERATOR.T_c_out )/PARAM.T0;
    y(6) = abs( x(2) - RECUPERATOR.T_h_out )/PARAM.T0;
    
    
    % combustor 3
    COMBUST3.FLUID = PARAM.FLUID;
    COMBUST3.T_in = RECUPERATOR.T_c_out;
    COMBUST3.p_in = RECUPERATOR.p_c_out;
    COMBUST3.p_out = COMBUST3.p_in;
    COMBUST3.mdot = RECUPERATOR.mdot_c;
    COMBUST3.Q = PARAM.Q_NG(3)/6.3*6.5;   % additional combustion mass flow rate (from Reference)
    COMBUST3.NTU = PARAM.NTU_combus3;
    [ COMBUST3.T_out, COMBUST3.X_dest ] = fcn_combustor( COMBUST3, PARAM );
    
    % RECUP turbine
    TURBINE_RECUP.FLUID = PARAM.FLUID;
    TURBINE_RECUP.eta_is = PARAM.eta_turbine;
    TURBINE_RECUP.p_in = COMBUST3.p_out;
    TURBINE_RECUP.p_out = PARAM.p0;
    TURBINE_RECUP.T_in = COMBUST3.T_out;
    TURBINE_RECUP.mdot = COMBUST3.mdot;
    [TURBINE_RECUP.W, TURBINE_RECUP.T_out, TURBINE_RECUP.X_dest] = fcn_turbine( TURBINE_RECUP, PARAM );
    
    y(8) = abs( RECUPERATOR.T_h_in - TURBINE_RECUP.T_out )/PARAM.T0;
    
    PLANT.W_t = TURBINE_HP.W + TURBINE_LP.W + TURBINE_RECUP.W;
    PLANT.W_net = TURBINE_HP.W + TURBINE_LP.W + TURBINE_RECUP.W - COMPRESSOR.W - PUMP.W ;
    PLANT.W_in = PARAM.e_liq*PUMP.mdot + 0.6*PARAM.Q_NG(1)/6.3*12.8;
    PLANT.eta_rt = PLANT.W_net/PLANT.W_in;

    s1a = py.CoolProp.CoolProp.PropsSI('Smass','P',PARAM.p0,'T',PARAM.T0,PREHEATER.FLUID_h);
    s2a = py.CoolProp.CoolProp.PropsSI('Smass','P',PREHEATER.p_h_out,'T',PREHEATER.T_h_out,PREHEATER.FLUID_h);
    s3a = py.CoolProp.CoolProp.PropsSI('Smass','P',COMPRESSOR.p_out,'T',COMPRESSOR.T_out,COMPRESSOR.FLUID);
    s4a = py.CoolProp.CoolProp.PropsSI('Smass','P',RECUPERATOR.p_h_out,'T',RECUPERATOR.T_h_out,RECUPERATOR.FLUID_h);
    s5a = py.CoolProp.CoolProp.PropsSI('Smass','P',COMBUST3.p_out,'T',COMBUST3.T_out,COMBUST3.FLUID);
    s6a = py.CoolProp.CoolProp.PropsSI('Smass','P',TURBINE_RECUP.p_out,'T',TURBINE_RECUP.T_out,TURBINE_RECUP.FLUID);
    s7a = py.CoolProp.CoolProp.PropsSI('Smass','P',RECUPERATOR.p_c_out,'T',RECUPERATOR.T_c_out,RECUPERATOR.FLUID_c);
    
    data_Ts_RAANG = [...
        s1a,PARAM.T0;
        s2a,PREHEATER.T_h_out;
        s3a,COMPRESSOR.T_out;
        s4a,RECUPERATOR.T_h_out;
        s5a,COMBUST3.T_out;
        s6a,TURBINE_RECUP.T_out;
        s7a,RECUPERATOR.T_c_out
        ];
    
catch
    
    y = 1000*ones(length(T),1);
    
end

