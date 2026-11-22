% driver simulate desired liquid air power plants (LAPPs). For more info.,
% please refer to https://smyng91.github.io/LiquidAirPlant/. If you
% use/reference the code, please cite the following paper:
%
% S. Yang, Solar-driven liquid air power plant modeling, design space exploration, and multi-objective optimization, Energy, 2022
% https://doi.org/10.1016/j.energy.2022.123324
%
% Also check out driver_verification for other configurations besides
% the solar-driven LAPP.
%
% Dependency: CoolProp for MATLAB
%
clear; clc; close all
addpath(fullfile(fileparts(mfilename('fullpath')), 'components'));
options = optimoptions('fsolve', 'Display', 'none', 'UseParallel', false, ...
    'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-12);

%%
% read input
plant_input;

% AAS config
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = solve_model(func, x0, options, 'AAS');
[~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
AAS
