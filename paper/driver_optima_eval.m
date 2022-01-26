clear; clc; close all

addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);

%% run
plant_input;



balanced = [...
    0.548587299087344;  
    1.076213198091668;  
    0.122806279847624;  
    4.172148090860937;  
    0.155894182073230;  
    3.739407410183461;  
    2.249099165017271;  
    4.327127891503268;  
    1.685345533978621 ];
eta_rt = [    0.5707    6.9120    0.1068    5.1406    0.5214    2.9179    2.8686    2.7397    3.4738    0.2071    0.5001 ];
eta_ex = [ 0.2266    0.3312    0.1445    9.5124    0.4933    2.2631    1.7217    3.4825    4.5327    0.1316    0.6082 ];

x =balanced;

PLANT.mdot_la = x(1);    
PLANT.mdot_preheater = x(2);  
PLANT.alpha_p_int = x(3);  
PTC.mdot = x(4);  
SHX.alpha1 = x(5);  
NTU.preheat =x(6);  
NTU.recup = x(7);  
NTU.shx1 = x(8);  
NTU.shx2 = x(9);  

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
AAS.W_net 
TURBINE_HP.W + TURBINE_LP.W

% AAS.COP
% plot_Ts;
