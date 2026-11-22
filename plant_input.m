% Input file of the ambient air/solar (AAS) liquid air power plant.
% Run from driver.m (the components folder must be on the MATLAB path).

% ambient
AMBIENT.T0 = 298;
AMBIENT.p0 = 0.1e6;
AMBIENT.u0 = 0.2;
AMBIENT.FLUID = 'AIR';
AMBIENT.T_sun = 5770;
AMBIENT.Gb = 800;
[AMBIENT.h0, AMBIENT.s0] = fcn_props(AMBIENT.FLUID, 'PT', AMBIENT.p0, AMBIENT.T0, 'H', 'S');

% liquid air storage
PLANT.FLUID = 'AIR';
PLANT.e_liq = 0.42*1e3*3600;
PLANT.T_store = 79;
PLANT.p_store = 0.2e6;

% plant
PLANT.mdot_la = 0.23;
PLANT.p_pump = 20e6;
PLANT.mdot_preheater = 0.33;
% PLANT.p_int = AMBIENT.p0 + (PLANT.p_pump-AMBIENT.p0)*.01;
PLANT.p_int = PLANT.p_pump/6.72;
[PLANT.h_la_in, PLANT.s_la_in] = fcn_props(PLANT.FLUID, 'PT', PLANT.p_store, PLANT.T_store, 'H', 'S');

% efficiency
PLANT.eta_pump = 0.75;
PLANT.eta_turbine = 0.75;
PLANT.eta_comp = 0.75;

% PTC
PTC.mdot = 9.51;
PTC.FLUID = 'INCOMP::S800';
PTC.eta_opt = 0.726;
% HTF loop inlet pressure. It must keep the oil liquid after the receiver
% pressure drop (~0.4 MPa here): S800 boils above 0.36 MPa at 550 K and
% 1.37 MPa at its 671 K limit, so 0.2 MPa is infeasible.
PTC.p_in = 2e6;
PTC.d_c = 0.115;
PTC.d_r = 0.07;
PTC.d_ri = PTC.d_r/1.05;    % receiver inner diameter
PTC.k_r = 280;              % receiver wall thermal conductivity [W/m/K]
PTC.roughness = 1e-3;       % receiver inner wall absolute roughness [m]
PTC.e_c = 0.86;
PTC.e_r = 0.14;
PTC.sigma = 5.67e-8;
PTC.rim = 70;
PTC.L = 122.7;
PTC.L_f = 1.84;
PTC.A_ap = 4*PTC.L_f*tand(PTC.rim/2)*PTC.L;
PTC.A_r = pi()*PTC.d_r*PTC.L;
PTC.A_c = pi()*PTC.d_c*PTC.L;
PTC.Cr = PTC.A_ap/PTC.A_r;

% SHX
SHX.alpha1 = 0.4932;
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;

% NTUs
NTU.total = 12;
NTU.preheat = 2.26;
NTU.recup = 1.72;
NTU.shx1 = 3.48;
NTU.shx2 = 4.54;
