clear; clc; close all
addpath('./components/');

% import design & operational parameters
plant_input;

%% sweep for eta_rt
Np = 20;
mdot_la = linspace(0.3,2.5,Np);
p_r = 7;
p2 = 20e6;
p0 = 0.1e6;

NTU.preheat = 3;
NTU.recup = 3;
NTU.shx1 = 3;
NTU.shx2 = 3;
PTC.mdot = 10;
SHX.alpha1 = 0.4;
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;
PLANT.mdot_preheater = 5;

idx = 1;
for i = 1:length(mdot_la)
    PLANT.mdot_la = mdot_la(i);
%     for j = 1:length(p_r)
        PLANT.p_int = p2/p_r;
        [f,g,h] = problem( AMBIENT, PLANT, NTU, PTC, SHX );
        solution(idx,:) = [ PLANT.mdot_la, f,g,h];
        idx=idx+1;
%     end
end
save('paper\solution_sweep_effectiveness.mat','solution');


%%

function [f,g,h] = problem( AMBIENT, PLANT, NTU, PTC, SHX )

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
    g = NaN*ones(1,9);
    h = NaN*ones(1,4);
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
    g(9) = RECUPERATOR.T_h_out;
    h(1) = PREHEATER.eff;
    h(2) = RECUPERATOR.eff;
    h(3) = SHX1.eff;
    h(4) = SHX2.eff;
end

end
