clear; clc; close all
addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);


%% verification
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