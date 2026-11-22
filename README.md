# LiquidAirPlant

[![View LiquidAirPlant on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/106080-liquidairplant)

## Overview 
LiquidAirPlant is a MATLAB code allowing its users to model various liquid air power plant configurations driven by natural gas combustors, parabolic trough solar collectors, and/or ambient air. 

## Citation
Please cite the following paper if you reference/use the code
> S. Yang, _Solar-driven liquid air power plant modeling, design space exploration, and multi-objective optimization_, Energy, 2022 [https://doi.org/10.1016/j.energy.2022.123324](https://doi.org/10.1016/j.energy.2022.123324)

* **Please note that the PTC rim angle in the paper should be 70 degrees, not 50 degrees. Credit to Kevin Lu at OSU for identifying/reporting the discrepancy.**

* **Feel free to use/improve the code/model as needed and submit pull requests.**

## Dependencies
* [MATLAB 2021a or above](https://www.mathworks.com/products/matlab.html)
* [CoolProp for MATLAB](http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html) (through MATLAB's Python interface, i.e. `py.CoolProp.CoolProp`). You may need to install CoolProp from source for Python 3.8 or above. Read [THIS](http://www.coolprop.org/coolprop/wrappers/Python/index.html) for manual installation.

## Quickstart
After cloning the repository, run
```MATLAB
driver_verification
```
in MATLAB to simulate ambient air-driven, natural gas-driven, and recuperative natural gas-driven liquid air power plants (LAPPs) as described in 

[Antonelli, M., Barsali, S., Desideri, U., Giglioli, R., Paganucci, F. and Pasini, G., 2017. Liquid air energy storage: Potential and challenges of hybrid power plants. Applied energy, 194, pp.522-529.](https://doi.org/10.1016/j.apenergy.2016.11.091)

which were referenced as validation cases for the proposed modeling framework.

You can also run 
```MATLAB
driver
```
to simulate a solar-driven LAPP as described in my [paper](https://doi.org/10.1016/j.energy.2022.123324). 

## Model description
The input file is
```MATLAB
plant_input.m
```
I recommend using structures for clarity/simplicity.

You can find all component models under
```MATLAB
.\components\
```
defined as MATLAB functions. You may also add new component models as functions, e.g., 
```MATLAB
function [out1, out2,...,out] = fcn_componentname( struct1, struct2,...,struct )
```

If you want to model and simulate other plant configurations, you can follow the same system assembly convention as
```MATLAB
model_X.m
```
where X is AA, AANG, RAANG, or AAS. 

Conventions used by the models:
* `y = model_X(T, ...)` returns only the signed, normalized residuals that `fsolve` drives to zero; post-processing (performance, exergy, T-s state points) runs only when more outputs are requested, e.g. `[y, PLANT] = model_X(T, ...)`.
* T-s state points are returned in `PLANT.Ts = [s, T]` (they were previously stored in the globals `data_Ts_X`).
* All CoolProp calls go through `components/fcn_props.m`, which reuses one CoolProp `AbstractState` per fluid so that a single flash calculation returns every property needed for a state.
* Components catch only CoolProp failures (physically infeasible states) and return penalty values so the solver can recover; other errors such as a missing input field are raised.
* `solve_model.m` wraps `fsolve` and warns when no steady state was found, so a non-converged result is never reported silently.

## Revision notes (2026)
The code was refactored for speed and robustness, and the following errors were corrected. Results therefore differ from the original release (commit `42a227a`):
* **Solver:** residuals were `abs()` values, which are non-differentiable at the solution and degrade `fsolve`'s finite-difference Jacobian. They are now signed and the solver tolerances are tighter; in testing, the residuals of the verification cases dropped from ~1e-5 to ~1e-11 with about half as many CoolProp evaluations.
* **Heat exchangers:** the effectiveness-NTU model used a point specific heat at the mean temperature, which violated the energy balance by up to ~10% of the duty for cryogenic air. It now uses the mean specific heat `dh/dT` over each stream's temperature span, so `mdot_c*dh_c = mdot_h*dh_h = Q` holds at convergence. The third output of `fcn_hx` is the duty `Q`; the effectiveness is now the fourth output (the AAS model previously stored `Q` as `eff`).
* **Combustors:** the outlet temperature now follows from the enthalpy balance instead of a constant inlet specific heat over a ~700 K temperature rise. The fuel mass now joins the working fluid (`mdot_fuel` input, `mdot_out` output; combustion products are treated as air), and the exergy destruction no longer carries an extra factor `T0`.
* **Fuel:** `LHV_NG` was 55.5 MJ/kg, which is the *higher* heating value of methane; it is now the LHV, 50.0 MJ/kg, consistent with the reference. Each combustor gets its own fuel flow (`mdot_combus1/2/3`), and `PLANT.Q_fuel` reports the total fuel heat input.
* **PTC:** the Haaland friction factor used `log` instead of `log10`; the convective coefficient, relative roughness and pressure drop divided by `1.05*d_r` instead of the inner diameter `d_r/1.05`; and a laminar friction branch was added. With the corrected (≈6× larger) pressure drop, the original loop pressure of 0.2 MPa is infeasible because S800 would boil in the receiver (even with the old pressure drop, the outlet pressure was below the S800 vapor pressure at the operating temperature). `PTC.p_in` is therefore 2 MPa, which keeps S800 liquid up to its 671 K limit. `PTC.d_ri`, `PTC.k_r` and `PTC.roughness` are now explicit inputs with the previous values.
* **AAS outputs:** `eta_II` (not the non-existent `eta_ex`) is set to NaN for wet turbine expansions, which are now detected from the actual (p, h) exhaust state. Failed component evaluations return NaN instead of zero power or exergy destruction.
* **Hard-coded values** (HP turbine outlet pressure, recuperative compressor pressure and fuel flow, liquid air specific energy in the AA model) moved to the input parameters with unchanged values.
* **Verification inputs** now follow the reference's stated conditions: fuel split 3.8 + 2.5 kg/h, complete combustion (`NTU_combus = Inf` instead of the former fitting factor 2.4), and liquid air warmed to 20 °C in the first ambient-air exchanger (`NTU_pre = 10.92`; the AA reheater is sized the same way).

## Verification against HYSYS
`driver_verification` prints a comparison with the HYSYS results in Table 1 of Antonelli et al., [Energy Procedia 88 (2016) 693-697](https://doi.org/10.1016/j.egypro.2016.06.100) (same authors and plant model as the 2017 paper, open access):

| Case | Quantity | HYSYS | This model | Difference |
|---|---|---|---|---|
| AA | pump power [kW] | 1.67 | 1.645 | -1.5% |
| AA | HP / LP turbine power [kW] | 4.61 / 7.35 | 4.698 / 7.449 | +1.9% / +1.4% |
| AA | net power [kW] | 10.29 | 10.502 | +2.1% |
| AANG | HP / LP turbine power [kW] | 32.73 / 50.35 | 32.832 / 50.230 | +0.3% / -0.2% |
| AANG | net power [kW] | 81.41 | 81.417 | matched (calibrated) |

* The reference does not state the size of the AANG exhaust regenerator. `NTU_reg = 13.8` (about a 60 K approach) reproduces its net power; the turbine split and the AA case are not calibrated. The resulting turbine inlet temperatures (1725 K and 1749 K) agree with those implied by the HYSYS turbine powers (about 1720 K and 1750 K). With `NTU_reg = 5` the AANG net power is 4.8% lower.
* The remaining 1-2% component-level differences are within the spread of property models: HYSYS used Peng-Robinson for an N2/O2 mixture, while this code uses CoolProp's reference equation of state for air. Treating the combustion products as air changes the turbine work by about 1% (checked against complete-combustion products with Cantera).
* No public reference values are available for RAANG and AAS.
* The MATLAB models were also cross-checked against an independent Python implementation, which agrees to better than 1e-9 in relative terms for all four configurations.

Round-trip efficiency reported by the models (equivalent definition with `eta_cc = 0.55`), original release → revised: AA 0.1005 → 0.1050, AANG 0.5300 → 0.5496, RAANG 0.6597 → 0.6718. The changes come from the corrections above and from the verification inputs that now follow the reference. For AANG the reference's own powers give 81.41/(100 + 0.55*87.5) = 0.5497.
