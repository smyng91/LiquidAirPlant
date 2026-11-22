% driver_verification verifies the liquid air power plant (LAPP) model
% proposed in
%
% S. Yang, Solar-driven liquid air power plant modeling, design space exploration, and multi-objective optimization, Energy, 2022
% https://doi.org/10.1016/j.energy.2022.123324
%
% against the HYSYS model published in
%
% Antonelli, M., Barsali, S., Desideri, U., Giglioli, R., Paganucci, F. and
% Pasini, G., 2017. Liquid air energy storage: Potential and challenges of
% hybrid power plants. Applied energy, 194, pp.522-529.
%
% The cited study evaluates several LAPP configurations and the proposed
% model verifies against three cases, namely ambient air-driven (AA),
% ambient air/natural gas-driven (AANG), and recuperative ambient
% air/natural gas-driven (RAANG).
%
% Dependency: CoolProp for MATLAB
% More info.: https://smyng91.github.io/LiquidAirPlant/
%
clear; clc; close all
addpath(fullfile(fileparts(mfilename('fullpath')), 'components'));
options = optimoptions('fsolve', 'Display', 'none', 'UseParallel', false, ...
    'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-12);

%% LAPP verification parameters
% Antonelli, M., Barsali, S., Desideri, U., Giglioli, R., Paganucci, F. and
% Pasini, G., 2017. Liquid air energy storage: Potential and challenges of
% hybrid power plants. Applied energy, 194, pp.522-529.

% ambient
PARAM.T0 = 298;
PARAM.p0 = 0.1e6;

% liquid air storage
PARAM.FLUID = 'AIR';
PARAM.e_liq = 0.5*1e3*3600;
PARAM.T_store = 78;
PARAM.p_store = 0.2e6;

% plant
PARAM.mdot_la = 200/3600;
PARAM.p_pump = 20e6;
PARAM.p_int = 3e6;          % HP turbine outlet pressure
PARAM.p_comp = 1.5e6;       % recuperative cycle compressor outlet pressure (RAANG)
PARAM.mdot_preheater = 0.12;
PARAM.mdot_regen = 1;
PARAM.LHV_NG = 55.5e6; % J/kg CH4 LHV
PARAM.mdot_combus = 6.3/3600;
PARAM.Q_NG =  PARAM.mdot_combus*PARAM.LHV_NG;
PARAM.mdot_combus3 = 6.5/3600;  % additional fuel of the recuperative cycle (RAANG)
PARAM.Q_NG3 = PARAM.mdot_combus3*PARAM.LHV_NG;

% NTUs
PARAM.NTU_pre = 5;
PARAM.NTU_reg = 5;
PARAM.NTU_rec = 5;
PARAM.NTU_combus1 = 2.4;
PARAM.NTU_combus2 = PARAM.NTU_combus1;
PARAM.NTU_combus3 = PARAM.NTU_combus1;

% efficiency
PARAM.eta_pump = 0.75;
PARAM.eta_turbine = 0.75;
PARAM.eta_comp = 0.75;
PARAM.eta_cc = 0.55;

%% Solve each configuration model
% AA config
func = @(T) model_AA( T, PARAM );
x0 = ones(4,1);
T = solve_model(func, x0, options, 'AA');
[~, AA] = model_AA( T, PARAM );
AA

% AANG config
func = @(T) model_AANG( T, PARAM );
x0 = ones(5,1);
T = solve_model(func, x0, options, 'AANG');
[~, AANG] = model_AANG( T, PARAM );
AANG

% RAANG config
func = @(T) model_RAANG( T, PARAM );
x0 = ones(8,1);
T = solve_model(func, x0, options, 'RAANG');
[~, RAANG] = model_RAANG( T, PARAM );
RAANG
