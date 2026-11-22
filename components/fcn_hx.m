function [y, X_dest, Q, effectiveness] = fcn_hx( HX, AMBIENT )
%FCN_HX Counter-flow heat exchanger, effectiveness-NTU method.
%
%   [y, X_dest, Q, effectiveness] = FCN_HX( HX, AMBIENT )
%
%   HX fields: FLUID_c, FLUID_h, mdot_c, mdot_h [kg/s], T_c_in, T_h_in [K],
%              p_c_in, p_h_in [Pa], NTU [-], and the solver's current guesses
%              of the outlet temperatures T_c_out, T_h_out [K]
%   y              computed outlet temperatures [T_c_out, T_h_out] [K]
%   X_dest         exergy destruction [W]
%   Q              heat duty [W]
%   effectiveness  [-]
%
%   The heat capacity rate of each stream uses its mean specific heat over
%   the actual temperature span, cp = (h_out - h_in)/(T_out - T_in), taken
%   at the guessed outlet temperatures. At the converged solution this makes
%   mdot_c*(h_c_out - h_c_in) = mdot_h*(h_h_in - h_h_out) = Q exactly, which
%   a point value of cp does not (the error reaches ~10% for cryogenic air
%   near its pseudo-critical region). Pressure losses are neglected.

% CoolProp failures (physically infeasible states) are caught so that the
% solver can recover; any other error is a genuine bug and is rethrown.
try
    [h_c_in, s_c_in] = fcn_props( HX.FLUID_c, 'PT', HX.p_c_in, HX.T_c_in, 'H', 'S' );
    [h_c_out, s_c_out] = fcn_props( HX.FLUID_c, 'PT', HX.p_c_in, HX.T_c_out, 'H', 'S' );
    [h_h_in, s_h_in] = fcn_props( HX.FLUID_h, 'PT', HX.p_h_in, HX.T_h_in, 'H', 'S' );
    [h_h_out, s_h_out] = fcn_props( HX.FLUID_h, 'PT', HX.p_h_in, HX.T_h_out, 'H', 'S' );

    Cc = HX.mdot_c*mean_cp( HX.FLUID_c, HX.p_c_in, HX.T_c_in, HX.T_c_out, h_c_in, h_c_out );
    Ch = HX.mdot_h*mean_cp( HX.FLUID_h, HX.p_h_in, HX.T_h_in, HX.T_h_out, h_h_in, h_h_out );
    Cmin = min(Cc, Ch);
    Cr = Cmin/max(Cc, Ch);

    NTU = HX.NTU;
    if abs(Cr - 1) < 1e-4
        effectiveness = NTU/(1 + NTU);
    elseif Cr < 1e-4
        effectiveness = 1 - exp(-NTU);
    else
        effectiveness = (1 - exp(-NTU*(1 - Cr)))/(1 - Cr*exp(-NTU*(1 - Cr)));
    end

    Q = effectiveness*Cmin*(HX.T_h_in - HX.T_c_in);
    y = [HX.T_c_in + Q/Cc, HX.T_h_in - Q/Ch];

    % exergy destruction, with the outlet states at the solver's outlet
    % temperatures (identical to the computed ones at convergence)
    X_dest = HX.mdot_c*( h_c_in - h_c_out - AMBIENT.T0*(s_c_in - s_c_out) ) + ...
        HX.mdot_h*( h_h_in - h_h_out - AMBIENT.T0*(s_h_in - s_h_out) );
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    y = [10000, 10000];
    X_dest = NaN;
    Q = NaN;
    effectiveness = NaN;
end

end

function cp = mean_cp( fluid, p, T_in, T_out, h_in, h_out )
% mean specific heat over [T_in, T_out]; point value if the span vanishes
if abs(T_out - T_in) > 1e-6*T_in
    cp = (h_out - h_in)/(T_out - T_in);
else
    cp = fcn_props( fluid, 'PT', p, T_in, 'CP' );
end
end
