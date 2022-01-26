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

solution = zeros(Np, 11);
eta_rt = zeros(Np,1); eta_II = zeros(Np,1); 
parfor i = 1:Np
    x = [ mdot_la(i), mdot_pre(i), alpha_p_int(i), mdot_ptc(i), alpha1(i),  NTU_pre(i), NTU_rec(i), NTU_shx1(i), NTU_shx2(i) ];
    [f,~] = fcn_evaluation( x, AMBIENT, PLANT, NTU, PTC, SHX );
    eta_rt(i) = -f(1);
    eta_II(i) = -f(2);
    sprintf('%i, %2.2f, %2.2f',i,-f(1),-f(2))
end
solution(:,1:9) = [mdot_la, mdot_pre, alpha_p_int, mdot_ptc, alpha1, NTU_pre, NTU_rec, NTU_shx1, NTU_shx2];
solution(:,10:11) = [eta_rt, eta_II];
save('paper\solution_sensitivity.mat','solution')