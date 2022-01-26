function [y, X_dest, Q] = fcn_hx( HX, AMBIENT )

% HX model based on effectiveness-NTU

% negligible pressure loss
p_c_out = HX.p_c_in;
p_h_out = HX.p_h_in;
T_h_avg = (HX.T_h_in + HX.T_h_out)/2;
T_c_avg = (HX.T_c_in + HX.T_c_out)/2;

% error/exception handling is needed to avoid divergence in physically
% infeasible design spaces, where CoolProp will output error and stop the
% code.
try
    
    cp_h = py.CoolProp.CoolProp.PropsSI('CPMASS','P',HX.p_h_in,'T',T_h_avg,HX.FLUID_h);
    cp_c = py.CoolProp.CoolProp.PropsSI('CPMASS','P',HX.p_c_in,'T',T_c_avg,HX.FLUID_c);
    Cc = HX.mdot_c*cp_c;
    Ch = HX.mdot_h*cp_h;
    Cmin = min(Cc,Ch);
    Cmax = max(Cc,Ch);
    Cr = Cmin/Cmax;
    
    if abs(Cr-1) < 1e-4
        effectiveness = HX.NTU/(1+HX.NTU);
    elseif Cr < 1e-4
        effectiveness = 1-exp(-HX.NTU);
    else
        effectiveness = (1-exp(-HX.NTU*(1-Cr)))/(1-Cr*exp(-HX.NTU*(1-Cr)));
    end
    
    dTmax = HX.T_h_in-HX.T_c_in;
    T_c_out = HX.T_c_in + Cmin/Cc*effectiveness*dTmax;
    T_h_out = HX.T_h_in - Cmin/Ch*effectiveness*dTmax;
    y = [T_c_out, T_h_out];
    
    h_h_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',HX.p_h_in,'T',HX.T_h_in,HX.FLUID_h);
    h_h_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',p_h_out,'T',T_h_out,HX.FLUID_h);
    h_c_in = py.CoolProp.CoolProp.PropsSI('Hmass','P',HX.p_c_in,'T',HX.T_c_in,HX.FLUID_c);
    h_c_out = py.CoolProp.CoolProp.PropsSI('Hmass','P',p_c_out,'T',T_c_out,HX.FLUID_c);
    
    s_h_in = py.CoolProp.CoolProp.PropsSI('Smass','P',HX.p_h_in,'T',HX.T_h_in,HX.FLUID_h);
    s_h_out = py.CoolProp.CoolProp.PropsSI('Smass','P',p_h_out,'T',T_h_out,HX.FLUID_h);
    s_c_in = py.CoolProp.CoolProp.PropsSI('Smass','P',HX.p_c_in,'T',HX.T_c_in,HX.FLUID_c);
    s_c_out = py.CoolProp.CoolProp.PropsSI('Smass','P',p_c_out,'T',T_c_out,HX.FLUID_c);
    X_dest = HX.mdot_c*(h_c_in-h_c_out-AMBIENT.T0*(s_c_in-s_c_out)) + ...
        HX.mdot_h*(h_h_in-h_h_out-AMBIENT.T0*(s_h_in-s_h_out));
    
    Q = effectiveness*Cmin*(HX.T_h_in-HX.T_c_in);
    
catch
    
    y(1) = 10000;
    y(2) = 10000;
    X_dest = 10000;
    effectiveness = 0;
    Q = 0;
    
end

end

