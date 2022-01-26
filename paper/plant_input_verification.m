% ambient
PARAM.T0 = 298;
PARAM.p0 = 0.1e6;

% liquid air storage
PARAM.FLUID = 'AIR';
PARAM.e_liq = 0.5*1e3*3600;
PARAM.T_store = 78;
PARAM.p_store = 0.2e6;

% plant
PARAM.mdot_la = 200/3600;
PARAM.p_pump = 20e6;
PARAM.mdot_preheater = 0.09;
PARAM.mdot_regen = 1;
PARAM.LHV_NG = 55.5e6; % J/kg CH4 LHV
PARAM.mdot_combus = 6.3/3600;
PARAM.Q_NG =  PARAM.mdot_combus*PARAM.LHV_NG*ones(3,1);

% NTUs
PARAM.NTU_pre = 5;
PARAM.NTU_reg = 5;
PARAM.NTU_rec = 5;
PARAM.NTU_combus1 = 2.5;
PARAM.NTU_combus2 = PARAM.NTU_combus1;
PARAM.NTU_combus3 = PARAM.NTU_combus1;

% efficiency
PARAM.eta_pump = 0.75;
PARAM.eta_turbine = 0.75;
PARAM.eta_comp = 0.75;

