function [dy, PERFORMANCE, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC ] = model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX )
%
% Liquid air plant model:
% Configuration - Ambient Air Solar
%

PREHEATER.T_c_out = T(1)*AMBIENT.T0;
PREHEATER.T_h_out = T(2)*AMBIENT.T0;
RECUPERATOR.T_c_out = T(3)*AMBIENT.T0;
RECUPERATOR.T_h_out = T(4)*AMBIENT.T0;
SHX1.T_c_out = T(5)*AMBIENT.T0;
SHX1.T_h_out = T(6)*AMBIENT.T0;
SHX2.T_c_out = T(7)*AMBIENT.T0;
SHX2.T_h_out = T(8)*AMBIENT.T0;
TURBINE_LP.T_out = T(9)*AMBIENT.T0;
PTC.T_r = T(10)*AMBIENT.T0;
PTC.T_c = T(11)*AMBIENT.T0;
PTC.T_out = T(12)*AMBIENT.T0;

% pump
PUMP.p_in = PLANT.p_store;
PUMP.T_in = PLANT.T_store;
PUMP.mdot = PLANT.mdot_la;
PUMP.eta_is = PLANT.eta_pump;
PUMP.p_out = PLANT.p_pump;
PUMP.FLUID = PLANT.FLUID;
[PUMP.W, PUMP.T_out, PUMP.X_dest] = fcn_pump( PUMP, AMBIENT );

% preheater
PREHEATER.FLUID_h = PLANT.FLUID;
PREHEATER.FLUID_c = PLANT.FLUID;
PREHEATER.mdot_c = PUMP.mdot;
PREHEATER.mdot_h = PLANT.mdot_preheater;
PREHEATER.T_h_in = AMBIENT.T0;
PREHEATER.p_h_in = AMBIENT.p0;
PREHEATER.T_c_in = PUMP.T_out;
PREHEATER.p_c_in = PUMP.p_out;
PREHEATER.p_c_out = PREHEATER.p_c_in;
PREHEATER.p_h_out = PREHEATER.p_h_in;
PREHEATER.NTU = NTU.preheat;
[x, PREHEATER.X_dest] = fcn_hx( PREHEATER, AMBIENT );
dy(1) = abs( x(1) - PREHEATER.T_c_out )/AMBIENT.T0;
dy(2) = abs( x(2) - PREHEATER.T_h_out )/AMBIENT.T0;

% recuperator
RECUPERATOR.FLUID_h = PLANT.FLUID;
RECUPERATOR.FLUID_c = PLANT.FLUID;
RECUPERATOR.mdot_c = PREHEATER.mdot_c;
RECUPERATOR.mdot_h = RECUPERATOR.mdot_c;
RECUPERATOR.T_c_in = PREHEATER.T_c_out;
RECUPERATOR.T_h_in = TURBINE_LP.T_out;
RECUPERATOR.p_c_in = PREHEATER.p_c_out;
RECUPERATOR.p_h_in = AMBIENT.p0;
RECUPERATOR.p_c_out = RECUPERATOR.p_c_in;
RECUPERATOR.p_h_out = RECUPERATOR.p_h_in;
RECUPERATOR.NTU = NTU.recup;
[x, RECUPERATOR.X_dest] = fcn_hx( RECUPERATOR, AMBIENT );
dy(3) = abs( x(1) - RECUPERATOR.T_c_out )/AMBIENT.T0;
dy(4) = abs( x(2) - RECUPERATOR.T_h_out )/AMBIENT.T0;

% PTC
SHX.T_mix = (SHX.mdot1*SHX1.T_h_out+SHX.mdot2*SHX2.T_h_out)/PTC.mdot;
PTC.T_in = SHX.T_mix;
[x, PTC.p_out, PTC.X_sun, PTC.X_dest, PTC.eta_c] = fcn_PTC( PTC, AMBIENT );
dy(10) = abs( x(1)-PTC.T_r )/AMBIENT.T0;
dy(11) = abs( x(2)-PTC.T_c )/AMBIENT.T0;
dy(12) = abs( x(3) - PTC.T_out )/AMBIENT.T0;

% solar HX 1
SHX1.FLUID_c = PLANT.FLUID;
SHX1.FLUID_h = PTC.FLUID;
SHX1.mdot_c = PUMP.mdot;
SHX1.mdot_h = SHX.mdot1;
SHX1.T_h_in = PTC.T_out;
SHX1.p_h_in = PTC.p_out;
SHX1.T_c_in = RECUPERATOR.T_c_out;
SHX1.p_c_in = RECUPERATOR.p_c_out;
SHX1.p_c_out = SHX1.p_c_in;
SHX1.p_h_out = SHX1.p_h_in;
SHX1.NTU = NTU.shx1;
[x, SHX1.X_dest] = fcn_hx( SHX1, AMBIENT );
dy(5) = abs( x(1) - SHX1.T_c_out )/AMBIENT.T0;
dy(6) = abs( x(2) - SHX1.T_h_out )/AMBIENT.T0;

% Turbine, HP
TURBINE_HP.FLUID = PLANT.FLUID;
TURBINE_HP.eta_is = PLANT.eta_turbine;
TURBINE_HP.p_in = SHX1.p_c_out;
TURBINE_HP.p_out = PLANT.p_int;
TURBINE_HP.T_in = SHX1.T_c_out;
TURBINE_HP.mdot = SHX1.mdot_c;
[TURBINE_HP.W, TURBINE_HP.T_out, TURBINE_HP.X_dest] = fcn_turbine( TURBINE_HP, AMBIENT );

% solar HX 2
SHX2.FLUID_c = PLANT.FLUID;
SHX2.FLUID_h = PTC.FLUID;
SHX2.mdot_c = TURBINE_HP.mdot;
SHX2.mdot_h = SHX.mdot2;
SHX2.T_h_in = PTC.T_out;
SHX2.p_h_in = PTC.p_out;
SHX2.T_c_in = TURBINE_HP.T_out;
SHX2.p_c_in = TURBINE_HP.p_out;
SHX2.p_c_out = SHX2.p_c_in;
SHX2.p_h_out = SHX2.p_h_in;
SHX2.NTU = NTU.shx2;
[x, SHX2.X_dest] = fcn_hx( SHX2, AMBIENT );
dy(7) = abs( x(1) - SHX2.T_c_out )/AMBIENT.T0;
dy(8) = abs( x(2) - SHX2.T_h_out )/AMBIENT.T0;

% Turbine, LP
TURBINE_LP.FLUID = PLANT.FLUID;
TURBINE_LP.eta_is = PLANT.eta_turbine;
TURBINE_LP.p_in = SHX2.p_c_out;
TURBINE_LP.p_out = AMBIENT.p0;
TURBINE_LP.T_in = SHX2.T_c_out;
TURBINE_LP.mdot = SHX2.mdot_c;
[TURBINE_LP.W, TURBINE_LP.T_out, TURBINE_LP.X_dest] = fcn_turbine( TURBINE_LP, AMBIENT );
dy(9) = abs( TURBINE_LP.T_out - RECUPERATOR.T_h_in )/AMBIENT.T0;

% PTC pump
PUMP_SHX.p_in = PTC.p_out;
PUMP_SHX.T_in = SHX.T_mix;
PUMP_SHX.mdot = PTC.mdot;
PUMP_SHX.eta_is = PLANT.eta_pump;
PUMP_SHX.p_out = PTC.p_in;
PUMP_SHX.FLUID = PTC.FLUID;
[PUMP_SHX.W, PUMP_SHX.T_out, PUMP_SHX.X_dest] = fcn_pump( PUMP_SHX, AMBIENT );

PERFORMANCE.W_net = TURBINE_HP.W + TURBINE_LP.W - PUMP.W - PUMP_SHX.W ;
PERFORMANCE.W_in = PLANT.e_liq*PUMP.mdot + AMBIENT.Gb*PTC.A_ap;

PERFORMANCE.X_in = PLANT.mdot_la*(PLANT.h_la_in-AMBIENT.h0 - AMBIENT.T0*(PLANT.s_la_in-AMBIENT.s0)) + PTC.X_sun;
PERFORMANCE.X_dest = [PUMP.X_dest, PREHEATER.X_dest, RECUPERATOR.X_dest,...
    TURBINE_HP.X_dest, TURBINE_LP.X_dest, SHX1.X_dest, SHX2.X_dest, PTC.X_dest, PUMP_SHX.X_dest ];

PERFORMANCE.eta_rt = PERFORMANCE.W_net / PERFORMANCE.W_in;
PERFORMANCE.eta_II = 1-sum(PERFORMANCE.X_dest)/PERFORMANCE.X_in;
PERFORMANCE.eta_c = PTC.eta_c;

try
    PREHEATER.h_h_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',AMBIENT.p0,'T',PREHEATER.T_h_out, AMBIENT.FLUID);
    RECUPERATOR.h_h_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',RECUPERATOR.p_h_in,'T',RECUPERATOR.T_h_out, PLANT.FLUID);
    PERFORMANCE.Qcool = PREHEATER.mdot_h*(AMBIENT.h0-PREHEATER.h_h_out);
    PERFORMANCE.COP = PERFORMANCE.Qcool/(PLANT.mdot_la*PLANT.e_liq+PUMP.W);
    phase_t1_out = py.CoolProp.CoolProp.PropsSI('Phase','T',TURBINE_HP.T_out,'P',TURBINE_HP.p_out,PLANT.FLUID);
    phase_t2_out = py.CoolProp.CoolProp.PropsSI('Phase','T',TURBINE_LP.T_out,'P',TURBINE_LP.p_out,PLANT.FLUID);
    if round(phase_t1_out) == 6 || round(phase_t2_out) == 6
        PERFORMANCE.eta_rt = NaN;
        PERFORMANCE.eta_ex = NaN;
        PERFORMANCE.Qcool = NaN;
        PERFORMANCE.COP = NaN;
    end
catch
    PERFORMANCE.eta_rt = NaN;
    PERFORMANCE.eta_ex = NaN;
    PERFORMANCE.Qcool = NaN;
    PERFORMANCE.COP = NaN;
end

end
