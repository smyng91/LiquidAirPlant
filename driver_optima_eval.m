clear; clc; close all

addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);

%% run
plant_input;

PLANT.mdot_la = 7.8286;    
PLANT.mdot_preheater = 1.4900;  
PLANT.alpha_p_int = 0.2072;  
PTC.mdot = 2.5836;  
SHX.alpha1 = 0.4239;  
NTU.preheat =7.1745;  
NTU.recup = 4.0566;  
NTU.shx1 = 0.4426;  
NTU.shx2 = 0.3262;  

PLANT.p_int = AMBIENT.p0 + (PLANT.p_pump-AMBIENT.p0)*PLANT.alpha_p_int;
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;

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
% AAS.COP
plot_Ts;
