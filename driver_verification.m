% driver_verification verifies the liquid air power plant (LAPP) model
% proposed in
% 
% paper
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

% read LAPP design parameters
plant_input_verification;

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