clear; clc; close all

addpath('./components/');
options = optimoptions('fsolve','Display','None',...
    'UseParallel',false);

%% run
% import design & operational parameters
plant_input;

% generate Sobol sequences
Np = 100000;
p = sobolset(9);

mdot_la = (1-0.01)*p(1:Np,1)+0.01;
mdot_preheater = (10-0.01)*p(1:Np,2)+0.01;
alpha_p_int = p(1:Np,3);
mdot_ptc = (20-0.1)*p(1:Np,4)+0.1;
shx_alpha1 = p(1:Np,5);
NTU_preheat = (5-0.1)*p(1:Np,6)+0.1;
NTU_recup = (5-0.1)*p(1:Np,7)+0.1;
NTU_shx1 = (5-0.1)*p(1:Np,8)+0.1;
NTU_shx2 = (5-0.1)*p(1:Np,9)+0.1;


sol = zeros(Np,12);
for i=1:Np
    
    PLANT.mdot_la = mdot_la(i);
    PLANT.mdot_preheater = mdot_preheater(i);
    PLANT.alpha_p_int = alpha_p_int(i);
    PLANT.p_int = AMBIENT.p0 + (PLANT.p_pump-AMBIENT.p0)*PLANT.alpha_p_int;
    PTC.mdot = mdot_ptc(i);
    SHX.alpha1 = shx_alpha1(i);
    SHX.mdot1 = PTC.mdot*SHX.alpha1;
    SHX.alpha2 = 1-SHX.alpha1;
    SHX.mdot2 = PTC.mdot*SHX.alpha2;
    NTU.preheat = NTU_preheat(i);
    NTU.recup = NTU_recup(i);
    NTU.shx1 = NTU_shx1(i);
    NTU.shx2 = NTU_shx2(i);
    
    % AAS config
    func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
    x0 = ones(12,1);
    [T, fval, exitflag, output] = fsolve(func, x0, options);
    [~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
        SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC] = ...
        model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
    
    if exitflag <= 0 || AAS.eta_II <= 0 || AAS.eta_rt <= 0 
        AAS.eta_II = NaN;
        AAS.eta_rt = NaN;
        Qcool = NaN;
        continue
    else
        PREHEATER.h_h_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',AMBIENT.p0,'T',PREHEATER.T_h_out, AMBIENT.FLUID);
        RECUPERATOR.h_h_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',RECUPERATOR.p_h_in,'T',RECUPERATOR.T_h_out, PLANT.FLUID);
        Qcool = PREHEATER.mdot_h*(AMBIENT.h0-PREHEATER.h_h_out) + ...
        RECUPERATOR.mdot_h*(AMBIENT.h0-RECUPERATOR.h_h_out);
    end
    sol(i,10) = AAS.eta_rt;
    sol(i,11) = AAS.eta_II;
    sol(i,12) = Qcool;
end


%%
sol(:,1) = mdot_la;
sol(:,2) = mdot_preheater;
sol(:,3) = alpha_p_int;
sol(:,4) = mdot_ptc;
sol(:,5) = shx_alpha1;
sol(:,6) = NTU_preheat;
sol(:,7) = NTU_recup;
sol(:,8) = NTU_shx1;
sol(:,9) = NTU_shx2;
save('solution_sensitivity.mat',sol);

% plot_Ts;
