clear; clc; close all

addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);

%% run
plant_input;

% AAS config
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
% AAS

AAS.eta_rt
AAS.eta_II
AAS.COP
% plot_Ts;
