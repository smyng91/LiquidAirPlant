% driver simulate desired liquid air power plants (LAPPs). For more info.,
% please refer to https://smyng91.github.io/LiquidAirPlant/. If you
% use/reference the code, please cite the following paper:
% 
% Energy
% 
% Also check out driver_verification for other configurations besides
% the solar-driven LAPP.
% 
% Dependency: CoolProp for MATLAB
% 
clear; clc; close all
addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);

%%
% read input
plant_input;

% AAS config
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
AAS
