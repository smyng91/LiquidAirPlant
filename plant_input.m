% ambient
AMBIENT.T0 = 298;
AMBIENT.p0 = 0.1e6;
AMBIENT.u0 = 0.2;
AMBIENT.FLUID = 'AIR';
AMBIENT.T_sun = 5770;
AMBIENT.Gb = 800;
AMBIENT.h0 = py.CoolProp.CoolProp.PropsSI('Hmass','P',AMBIENT.p0,'T',AMBIENT.T0, AMBIENT.FLUID);
AMBIENT.s0 = py.CoolProp.CoolProp.PropsSI('Smass','P',AMBIENT.p0,'T',AMBIENT.T0, AMBIENT.FLUID);

% liquid air storage
PLANT.FLUID = 'AIR';
PLANT.e_liq = 0.42*1e3*3600;
PLANT.T_store = 79;
PLANT.p_store = 0.2e6;

% plant
PLANT.mdot_la = 0.6;
PLANT.p_pump = 20e6;
PLANT.mdot_preheater = 5;
PLANT.p_int = AMBIENT.p0 + (PLANT.p_pump-AMBIENT.p0)*.01;
PLANT.h_la_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',PLANT.p_store,'T',PLANT.T_store, PLANT.FLUID);
PLANT.s_la_in = py.CoolProp.CoolProp.PropsSI('Smass','P',PLANT.p_store,'T',PLANT.T_store, PLANT.FLUID);

% efficiency
PLANT.eta_pump = 0.75;
PLANT.eta_turbine = 0.75;
PLANT.eta_comp = 0.75;

% PTC
PTC.mdot = 10;
PTC.FLUID = 'INCOMP::S800';
PTC.eta_opt = 0.726;
PTC.p_in = 0.2e6;
PTC.d_c = 0.115;
PTC.d_r = 0.07;
PTC.e_c = 0.86;
PTC.e_r = 0.14;
PTC.sigma = 5.67e-8;
PTC.rim = 50;
PTC.L = 122.7;
PTC.L_f = 1.84;
PTC.A_ap = 4*PTC.L_f*tand(PTC.rim/2)*PTC.L;
PTC.A_r = pi()*PTC.d_r*PTC.L; 
PTC.A_c = pi()*PTC.d_c*PTC.L; 
PTC.Cr = PTC.A_ap/PTC.A_r;

% SHX
SHX.alpha1 = 0.4;
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;

% NTUs
NTU.total = 12;
NTU.preheat = 3;
NTU.recup = 3;
NTU.shx1 = 3;
NTU.shx2 = 3;

    

