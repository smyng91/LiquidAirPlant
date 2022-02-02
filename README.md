# LiquidAirPlant

[![View LiquidAirPlant on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/106080-liquidairplant)

## Overview 
LiquidAirPlant is a MATLAB code allowing its users to model various liquid air power plant configurations driven by natural gas combustors, parabolic trough solar collectors, and/or ambient air. 

## Citation
Please cite the following paper if you reference/use the code
> S. Yang, _Solar-driven liquid air power plant modeling, design space exploration, and multi-objective optimization_, Energy, 2022 [https://doi.org/10.1016/j.energy.2022.123324](https://doi.org/10.1016/j.energy.2022.123324)

**Feel free to use/improve the code/model as needed and submit pull requests.**

## Dependencies
* [MATLAB 2021a or above](https://www.mathworks.com/products/matlab.html)
* [CoolProp for MATLAB](http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html)

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
where x is AA, AANG, RAANG, or AAS. 
