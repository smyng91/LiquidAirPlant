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
% Antonelli, M., Desideri, U., Giglioli, R., Paganucci, F. and Pasini, G.,
% 2016. Liquid air energy storage: a potential low emissions and efficient
% storage system. Energy Procedia, 88, pp.693-697 (open access; Table 1
% lists the HYSYS results of the AA and AANG cases used for comparison below).
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
% Antonelli et al. (2016, 2017), see above

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

% natural gas (methane), fuel split between the combustors as in the
% reference (6.3 kg/h = 3.8 + 2.5 kg/h)
PARAM.LHV_NG = 50.0e6;          % J/kg, CH4 lower heating value
PARAM.mdot_combus1 = 3.8/3600;  % combustor 1 (before the HP turbine)
PARAM.mdot_combus2 = 2.5/3600;  % combustor 2 (before the LP turbine)
PARAM.mdot_combus3 = 6.5/3600;  % combustor 3, recuperative cycle (RAANG)

% NTUs
PARAM.NTU_pre = 10.92;      % sized so the liquid air leaves at 20 C, as in the reference
PARAM.NTU_reg = 13.8;       % exhaust regenerator, not specified in the reference: sized to
                            % reproduce its AANG net power (approach ~60 K); the HP/LP
                            % turbine split then serves as an independent check
PARAM.NTU_rec = 5;
PARAM.NTU_combus1 = Inf;    % complete, adiabatic heat release (HYSYS conversion reactors)
PARAM.NTU_combus2 = PARAM.NTU_combus1;
PARAM.NTU_combus3 = PARAM.NTU_combus1;

% efficiency
PARAM.eta_pump = 0.75;
PARAM.eta_turbine = 0.75;
PARAM.eta_comp = 0.75;
PARAM.eta_cc = 0.55;

% AA: the ambient-air reheater after the HP turbine also delivers 20 C
% (assumed; the reference states 20 C for the first ambient-air exchanger)
PARAM_AA = PARAM;
PARAM_AA.NTU_reg = 3.25;

%% Solve each configuration model
% AA config
func = @(T) model_AA( T, PARAM_AA );
x0 = ones(4,1);
T = solve_model(func, x0, options, 'AA');
[~, AA, COMP_AA] = model_AA( T, PARAM_AA );
AA

% AANG config
func = @(T) model_AANG( T, PARAM );
x0 = ones(5,1);
T = solve_model(func, x0, options, 'AANG');
[~, AANG, COMP_AANG] = model_AANG( T, PARAM );
AANG

% RAANG config
func = @(T) model_RAANG( T, PARAM );
x0 = ones(8,1);
T = solve_model(func, x0, options, 'RAANG');
[~, RAANG, COMP_RAANG] = model_RAANG( T, PARAM );
RAANG

%% Comparison with HYSYS, Antonelli et al. (2016), Table 1
% columns: pump, HP turbine, LP turbine and net power [kW]
HYSYS.AA = [1.67, 4.61, 7.35, 10.29];
HYSYS.AANG = [1.67, 32.73, 50.35, 81.41];
MODEL.AA = [COMP_AA.PUMP.W, COMP_AA.TURBINE_HP.W, COMP_AA.TURBINE_LP.W, AA.W_net]/1e3;
MODEL.AANG = [COMP_AANG.PUMP.W, COMP_AANG.TURBINE_HP.W, COMP_AANG.TURBINE_LP.W, AANG.W_net]/1e3;
% round-trip efficiency as defined in the reference, P_net/(E_air*mdot + fuel LHV)
E_air = PARAM.e_liq*PARAM.mdot_la/1e3;
Q_fuel = (PARAM.mdot_combus1 + PARAM.mdot_combus2)*PARAM.LHV_NG/1e3;
HYSYS.AA(5) = HYSYS.AA(4)/E_air;                MODEL.AA(5) = MODEL.AA(4)/E_air;
HYSYS.AANG(5) = HYSYS.AANG(4)/(E_air + Q_fuel); MODEL.AANG(5) = MODEL.AANG(4)/(E_air + Q_fuel);

names = {'P_pump [kW]', 'P_HP [kW]', 'P_LP [kW]', 'P_net [kW]', 'eta_rt [-]'};
fprintf('\nComparison with HYSYS (Antonelli et al., Energy Procedia 88, 2016, Table 1)\n');
fprintf('%-6s %-12s %9s %9s %8s\n', 'case', 'quantity', 'HYSYS', 'model', 'diff');
cases = {'AA', 'AANG'};
for i = 1:numel(cases)
    for k = 1:numel(names)
        ref = HYSYS.(cases{i})(k); val = MODEL.(cases{i})(k);
        fprintf('%-6s %-12s %9.3f %9.3f %+7.1f%%\n', cases{i}, names{k}, ref, val, 100*(val/ref - 1));
    end
end
fprintf('AANG turbine inlet temperatures: HP %.0f K, LP %.0f K\n', ...
    COMP_AANG.TURBINE_HP.T_in, COMP_AANG.TURBINE_LP.T_in);
fprintf(['Note: AANG P_net is matched by sizing the regenerator (NTU_reg); the ', ...
    'other AANG rows and all AA rows are uncalibrated.\n']);
