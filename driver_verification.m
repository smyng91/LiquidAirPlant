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
addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);

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
PARAM.mdot_preheater = 0.12;
PARAM.mdot_regen = 1;
PARAM.LHV_NG = 55.5e6; % J/kg CH4 LHV
PARAM.mdot_combus = 6.3/3600;
PARAM.Q_NG =  PARAM.mdot_combus*PARAM.LHV_NG;

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
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AA] = model_AA( T, PARAM );
AA

% AANG config
func = @(T) model_AANG( T, PARAM );
x0 = ones(5,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AANG] = model_AANG( T, PARAM );
AANG

% RAANG config
func = @(T) model_RAANG( T, PARAM );
x0 = ones(8,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, RAANG] = model_RAANG( T, PARAM );
RAANG