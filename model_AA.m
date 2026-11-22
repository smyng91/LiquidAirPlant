function [y, PLANT, COMP] = model_AA( T, PARAM )
%MODEL_AA Liquid air power plant, ambient air-driven (AA) configuration.
%
%   y = MODEL_AA( T, PARAM ) returns the residuals for fsolve. T holds the
%   unknown stream temperatures normalized by PARAM.T0:
%     T(1:2)  preheater cold/hot outlet
%     T(3:4)  regenerator cold/hot outlet
%
%   [y, PLANT, COMP] = MODEL_AA( T, PARAM ) also returns the plant
%   performance (including the T-s state points in PLANT.Ts = [s, T]) and
%   the component structs. The post-processing is skipped when only the
%   residuals are requested, which keeps the solver iterations cheap.

T0 = PARAM.T0;
full = nargout > 1;
y = zeros(4,1);

% pump
PUMP.FLUID = PARAM.FLUID;
PUMP.p_in = PARAM.p_store;
PUMP.T_in = PARAM.T_store;
PUMP.mdot = PARAM.mdot_la;
PUMP.eta_is = PARAM.eta_pump;
PUMP.p_out = PARAM.p_pump;
[PUMP.W, PUMP.T_out, PUMP.X_dest] = fcn_pump( PUMP, PARAM );

% preheater (liquid air heated by ambient air)
PREHEATER.FLUID_h = PARAM.FLUID;
PREHEATER.FLUID_c = PARAM.FLUID;
PREHEATER.mdot_c = PUMP.mdot;
PREHEATER.mdot_h = PARAM.mdot_preheater;
PREHEATER.T_h_in = T0;
PREHEATER.T_c_in = PUMP.T_out;
PREHEATER.p_c_in = PUMP.p_out;
PREHEATER.p_h_in = PARAM.p0;
PREHEATER.p_c_out = PREHEATER.p_c_in;
PREHEATER.p_h_out = PREHEATER.p_h_in;
PREHEATER.NTU = PARAM.NTU_pre;
PREHEATER.T_c_out = T(1)*T0;
PREHEATER.T_h_out = T(2)*T0;
[x, PREHEATER.X_dest, PREHEATER.Q, PREHEATER.eff] = fcn_hx( PREHEATER, PARAM );
y(1) = ( x(1) - PREHEATER.T_c_out )/T0;
y(2) = ( x(2) - PREHEATER.T_h_out )/T0;

% HP turbine
TURBINE_HP.FLUID = PARAM.FLUID;
TURBINE_HP.eta_is = PARAM.eta_turbine;
TURBINE_HP.p_in = PREHEATER.p_c_out;
TURBINE_HP.p_out = PARAM.p_int;
TURBINE_HP.T_in = PREHEATER.T_c_out;
TURBINE_HP.mdot = PREHEATER.mdot_c;
[TURBINE_HP.W, TURBINE_HP.T_out, TURBINE_HP.X_dest] = fcn_turbine( TURBINE_HP, PARAM );

% regenerator (reheat by ambient air)
REGENERATOR.FLUID_h = PARAM.FLUID;
REGENERATOR.FLUID_c = PARAM.FLUID;
REGENERATOR.mdot_c = TURBINE_HP.mdot;
REGENERATOR.mdot_h = PARAM.mdot_regen;
REGENERATOR.T_h_in = T0;
REGENERATOR.T_c_in = TURBINE_HP.T_out;
REGENERATOR.p_c_in = TURBINE_HP.p_out;
REGENERATOR.p_h_in = PARAM.p0;
REGENERATOR.p_c_out = REGENERATOR.p_c_in;
REGENERATOR.p_h_out = REGENERATOR.p_h_in;
REGENERATOR.NTU = PARAM.NTU_reg;
REGENERATOR.T_c_out = T(3)*T0;
REGENERATOR.T_h_out = T(4)*T0;
[x, REGENERATOR.X_dest, REGENERATOR.Q, REGENERATOR.eff] = fcn_hx( REGENERATOR, PARAM );
y(3) = ( x(1) - REGENERATOR.T_c_out )/T0;
y(4) = ( x(2) - REGENERATOR.T_h_out )/T0;

if ~full, return; end

% LP turbine
TURBINE_LP.FLUID = PARAM.FLUID;
TURBINE_LP.eta_is = PARAM.eta_turbine;
TURBINE_LP.p_in = REGENERATOR.p_c_out;
TURBINE_LP.p_out = PARAM.p0;
TURBINE_LP.T_in = REGENERATOR.T_c_out;
TURBINE_LP.mdot = REGENERATOR.mdot_c;
[TURBINE_LP.W, TURBINE_LP.T_out, TURBINE_LP.X_dest] = fcn_turbine( TURBINE_LP, PARAM );

% performance
PLANT.W_t = TURBINE_HP.W + TURBINE_LP.W;
PLANT.W_net = PLANT.W_t - PUMP.W;
PLANT.W_in = PARAM.e_liq*PUMP.mdot;
PLANT.eta_rt = PLANT.W_net/PLANT.W_in;

% T-s diagram state points [s, T]
PLANT.Ts = fcn_Ts( PARAM.FLUID, ...
    [PUMP.p_in, PUMP.p_out, PREHEATER.p_c_out, TURBINE_HP.p_out, REGENERATOR.p_c_out, TURBINE_LP.p_out], ...
    [PUMP.T_in, PUMP.T_out, PREHEATER.T_c_out, TURBINE_HP.T_out, REGENERATOR.T_c_out, TURBINE_LP.T_out] );

COMP = struct('PUMP', PUMP, 'PREHEATER', PREHEATER, 'TURBINE_HP', TURBINE_HP, ...
    'REGENERATOR', REGENERATOR, 'TURBINE_LP', TURBINE_LP);

end
