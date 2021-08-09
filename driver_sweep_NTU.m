clear; clc; close all
addpath('./components/');

% import design & operational parameters
plant_input;

%% sweep for eta_rt
Np = 40;
mdot_la = 0.5;
p_r = 4;
p2 = 20e6;
p0 = 0.1e6;

NTU_preheat = linspace(0.01,5.99,Np);
NTU.recup = 3;
NTU.shx2 = 3;
PTC.mdot = 10;
SHX.alpha1 = 0.4;
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;
PLANT.mdot_preheater = 5;
PLANT.p_int = p2/p_r;
PLANT.mdot_la = mdot_la;

idx = 1;
for i = 1:length(NTU_preheat)
    NTU.preheat = NTU_preheat(i);
    NTU.shx1 = 6-NTU_preheat(i);
    [f,g] = problem( AMBIENT, PLANT, NTU, PTC, SHX );
    solution(idx,:) = [ NTU.preheat, NTU.shx1, f, g];
    idx=idx+1;
end
save('paper\solution_sweep_NTU.mat','solution');


%%

function [f,g] = problem( AMBIENT, PLANT, NTU, PTC, SHX )

% AAS config
options = optimoptions('fsolve','Display','None','UseParallel',false);
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC ] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );

if exitflag <= 0 || AAS.eta_II <= 0 || AAS.eta_rt <= 0
    f = NaN*ones(1,12);
    g = NaN*ones(1,8);
else
    % Objective functions F(X)
    f = [AAS.eta_rt, AAS.eta_II, AAS.COP, ...
        PUMP.X_dest, PUMP_SHX.X_dest, PREHEATER.X_dest, RECUPERATOR.X_dest, ...
        SHX1.X_dest, SHX2.X_dest, TURBINE_HP.X_dest, TURBINE_LP.X_dest, PTC.X_dest];
    g(1) = TURBINE_HP.T_in;
    g(2) = TURBINE_HP.T_out;
    g(3) = TURBINE_LP.T_out;
    g(4) = PREHEATER.T_c_out;
    g(5) = RECUPERATOR.T_c_out;
    g(6) = SHX.T_mix;
    g(7) = PTC.T_out;
    g(8) = TURBINE_LP.T_in;
end

end
