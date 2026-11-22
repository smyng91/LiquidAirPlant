function [dy, PERFORMANCE, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC ] = model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX )
%MODEL_AAS Liquid air power plant, ambient air/solar (AAS) configuration.
%
%   dy = MODEL_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX ) returns the residuals
%   for fsolve. T holds the unknown temperatures normalized by AMBIENT.T0:
%     T(1:2)    preheater cold/hot outlet
%     T(3:4)    recuperator cold/hot outlet
%     T(5:6)    solar HX 1 cold/hot outlet
%     T(7:8)    solar HX 2 cold/hot outlet
%     T(9)      LP turbine outlet (recuperator hot inlet)
%     T(10:12)  PTC receiver, glass cover and HTF outlet
%
%   [dy, PERFORMANCE, PUMP, ...] = MODEL_AAS( ... ) also returns the plant
%   performance and the component structs. The post-processing is skipped
%   when only the residuals are requested, which keeps the solver
%   iterations cheap.

T0 = AMBIENT.T0;
full = nargout > 1;
dy = zeros(12,1);

% pump
PUMP.p_in = PLANT.p_store;
PUMP.T_in = PLANT.T_store;
PUMP.mdot = PLANT.mdot_la;
PUMP.eta_is = PLANT.eta_pump;
PUMP.p_out = PLANT.p_pump;
PUMP.FLUID = PLANT.FLUID;
[PUMP.W, PUMP.T_out, PUMP.X_dest] = fcn_pump( PUMP, AMBIENT );

% preheater (liquid air cools ambient air)
PREHEATER.FLUID_h = PLANT.FLUID;
PREHEATER.FLUID_c = PLANT.FLUID;
PREHEATER.mdot_c = PUMP.mdot;
PREHEATER.mdot_h = PLANT.mdot_preheater;
PREHEATER.T_h_in = T0;
PREHEATER.p_h_in = AMBIENT.p0;
PREHEATER.T_c_in = PUMP.T_out;
PREHEATER.p_c_in = PUMP.p_out;
PREHEATER.p_c_out = PREHEATER.p_c_in;
PREHEATER.p_h_out = PREHEATER.p_h_in;
PREHEATER.NTU = NTU.preheat;
PREHEATER.T_c_out = T(1)*T0;
PREHEATER.T_h_out = T(2)*T0;
[x, PREHEATER.X_dest, PREHEATER.Q, PREHEATER.eff] = fcn_hx( PREHEATER, AMBIENT );
dy(1) = ( x(1) - PREHEATER.T_c_out )/T0;
dy(2) = ( x(2) - PREHEATER.T_h_out )/T0;

% recuperator (heated by the LP turbine exhaust)
TURBINE_LP.T_out = T(9)*T0;
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
RECUPERATOR.T_c_out = T(3)*T0;
RECUPERATOR.T_h_out = T(4)*T0;
[x, RECUPERATOR.X_dest, RECUPERATOR.Q, RECUPERATOR.eff] = fcn_hx( RECUPERATOR, AMBIENT );
dy(3) = ( x(1) - RECUPERATOR.T_c_out )/T0;
dy(4) = ( x(2) - RECUPERATOR.T_h_out )/T0;

% PTC, fed by the mixed HTF returning from both solar HXs
SHX1.T_h_out = T(6)*T0;
SHX2.T_h_out = T(8)*T0;
SHX.T_mix = (SHX.mdot1*SHX1.T_h_out + SHX.mdot2*SHX2.T_h_out)/PTC.mdot;
PTC.T_in = SHX.T_mix;
PTC.T_r = T(10)*T0;
PTC.T_c = T(11)*T0;
PTC.T_out = T(12)*T0;
[x, PTC.p_out, PTC.X_sun, PTC.X_dest, PTC.eta_c, PTC.Qu, PTC.dp] = fcn_PTC( PTC, AMBIENT );
dy(10) = ( x(1) - PTC.T_r )/T0;
dy(11) = ( x(2) - PTC.T_c )/T0;
dy(12) = ( x(3) - PTC.T_out )/T0;

% solar HX 1 (before the HP turbine)
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
SHX1.T_c_out = T(5)*T0;
[x, SHX1.X_dest, SHX1.Q, SHX1.eff] = fcn_hx( SHX1, AMBIENT );
dy(5) = ( x(1) - SHX1.T_c_out )/T0;
dy(6) = ( x(2) - SHX1.T_h_out )/T0;

% HP turbine
TURBINE_HP.FLUID = PLANT.FLUID;
TURBINE_HP.eta_is = PLANT.eta_turbine;
TURBINE_HP.p_in = SHX1.p_c_out;
TURBINE_HP.p_out = PLANT.p_int;
TURBINE_HP.T_in = SHX1.T_c_out;
TURBINE_HP.mdot = SHX1.mdot_c;
[TURBINE_HP.W, TURBINE_HP.T_out, TURBINE_HP.X_dest, TURBINE_HP.phase] = fcn_turbine( TURBINE_HP, AMBIENT );

% solar HX 2 (reheat before the LP turbine)
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
SHX2.T_c_out = T(7)*T0;
[x, SHX2.X_dest, SHX2.Q, SHX2.eff] = fcn_hx( SHX2, AMBIENT );
dy(7) = ( x(1) - SHX2.T_c_out )/T0;
dy(8) = ( x(2) - SHX2.T_h_out )/T0;

% LP turbine
TURBINE_LP.FLUID = PLANT.FLUID;
TURBINE_LP.eta_is = PLANT.eta_turbine;
TURBINE_LP.p_in = SHX2.p_c_out;
TURBINE_LP.p_out = AMBIENT.p0;
TURBINE_LP.T_in = SHX2.T_c_out;
TURBINE_LP.mdot = SHX2.mdot_c;
[TURBINE_LP.W, TURBINE_LP.T_out, TURBINE_LP.X_dest, TURBINE_LP.phase] = fcn_turbine( TURBINE_LP, AMBIENT );
dy(9) = ( TURBINE_LP.T_out - RECUPERATOR.T_h_in )/T0;

if ~full, return; end

% PTC pump (makes up the receiver pressure drop)
PUMP_SHX.p_in = PTC.p_out;
PUMP_SHX.T_in = SHX.T_mix;
PUMP_SHX.mdot = PTC.mdot;
PUMP_SHX.eta_is = PLANT.eta_pump;
PUMP_SHX.p_out = PTC.p_in;
PUMP_SHX.FLUID = PTC.FLUID;
[PUMP_SHX.W, PUMP_SHX.T_out, PUMP_SHX.X_dest] = fcn_pump( PUMP_SHX, AMBIENT );

% performance
PERFORMANCE.W_net = TURBINE_HP.W + TURBINE_LP.W - PUMP.W - PUMP_SHX.W;
PERFORMANCE.W_in = PLANT.e_liq*PUMP.mdot + AMBIENT.Gb*PTC.A_ap;
PERFORMANCE.X_in = PLANT.mdot_la*( PLANT.h_la_in - AMBIENT.h0 - T0*(PLANT.s_la_in - AMBIENT.s0) ) + PTC.X_sun;
PERFORMANCE.X_dest = [PUMP.X_dest, PREHEATER.X_dest, RECUPERATOR.X_dest, ...
    TURBINE_HP.X_dest, TURBINE_LP.X_dest, SHX1.X_dest, SHX2.X_dest, PTC.X_dest, PUMP_SHX.X_dest];
PERFORMANCE.eta_rt = PERFORMANCE.W_net/PERFORMANCE.W_in;
PERFORMANCE.eta_II = 1 - sum(PERFORMANCE.X_dest)/PERFORMANCE.X_in;
PERFORMANCE.eta_c = PTC.eta_c;

% cooling delivered by the preheater
try
    PREHEATER.h_h_out = fcn_props( AMBIENT.FLUID, 'PT', AMBIENT.p0, PREHEATER.T_h_out, 'H' );
    RECUPERATOR.h_h_out = fcn_props( PLANT.FLUID, 'PT', RECUPERATOR.p_h_in, RECUPERATOR.T_h_out, 'H' );
    PERFORMANCE.Qcool = PREHEATER.mdot_h*(AMBIENT.h0 - PREHEATER.h_h_out);
    PERFORMANCE.COP = PERFORMANCE.Qcool/(PLANT.mdot_la*PLANT.e_liq + PUMP.W);
    cooling_ok = true;
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    cooling_ok = false;
end

% reject designs with a wet (two-phase, CoolProp phase index 6) or failed
% turbine expansion, or an unevaluable preheater outlet
valid = cooling_ok && ...
    all(isfinite([TURBINE_HP.phase, TURBINE_LP.phase])) && ...
    TURBINE_HP.phase ~= 6 && TURBINE_LP.phase ~= 6;
if ~valid
    PERFORMANCE.eta_rt = NaN;
    PERFORMANCE.eta_II = NaN;
    PERFORMANCE.Qcool = NaN;
    PERFORMANCE.COP = NaN;
end

end
