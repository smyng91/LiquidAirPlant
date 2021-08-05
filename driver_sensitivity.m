clear; clc; close all
addpath('./components/');

% import design & operational parameters
plant_input;

Np = 100000;
p = sobolset(6);
mdot_la = (10-0.01)*p(1:Np,1)+0.01;
mdot_pre = (10-0.01)*p(1:Np,2)+0.01;
alpha_p_int = (0.99-0.01)*p(1:Np,3)+0.01;
mdot_ptc = (10-0.01)*p(1:Np,4)+0.01;
alpha1 = (0.99-0.01)*p(1:Np,5)+0.01;
r = rand(Np, 4);
r = r./sum(r,2);
NTU_tot = 12;
NTU_pre = r(:,1)*NTU_tot;
NTU_rec = r(:,2)*NTU_tot;
NTU_shx1 = r(:,3)*NTU_tot;
NTU_shx2 = r(:,4)*NTU_tot;

solution = zeros(Np, 12);
eta_rt = zeros(Np,1); eta_II = zeros(Np,1); COP = zeros(Np,1);
parfor i = 1:Np
    x = [ mdot_la(i), mdot_pre(i), alpha_p_int(i), mdot_ptc(i), alpha1(i),  NTU_pre(i), NTU_rec(i), NTU_shx1(i), NTU_shx2(i) ];
    f = problem( x, AMBIENT, PLANT, NTU, PTC, SHX );
    eta_rt(i) = f(1);
    eta_II(i) = f(2);
    COP(i)  = f(3);
end
solution(:,1:9) = [mdot_la, mdot_pre, alpha_p_int, mdot_ptc, alpha1, NTU_pre, NTU_rec, NTU_shx1, NTU_shx2];
solution(:,10:12) = [eta_rt, eta_II, COP];
 save('solution.mat','solution')

function f = problem( x, AMBIENT, PLANT, NTU, PTC, SHX )

PLANT.mdot_la = x(1);
PLANT.mdot_preheater = x(2);
PLANT.alpha_p_int = x(3);
PLANT.p_int = AMBIENT.p0 + (PLANT.p_pump-AMBIENT.p0)*PLANT.alpha_p_int;
PTC.mdot = x(4);
SHX.alpha1 = x(5);
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;
NTU.preheat = x(6);
NTU.recup = x(7);
NTU.shx1 = x(8);
NTU.shx2 = x(9);

% AAS config
options = optimoptions('fsolve','Display','None','UseParallel',false);
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AAS] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );

if exitflag <= 0 || AAS.eta_II <= 0 || AAS.eta_rt <= 0
    AAS.eta_II = NaN;
    AAS.eta_rt = NaN;
    AAS.COP = NaN;
end

% Objective functions F(X)
f(1) = AAS.eta_rt;
f(2) = AAS.eta_II;
f(3) = AAS.COP;

end
