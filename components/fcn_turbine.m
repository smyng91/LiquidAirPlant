function [W, T_out, X_dest, phase_out] = fcn_turbine( TURBINE, AMBIENT )
%FCN_TURBINE Adiabatic turbine with isentropic efficiency.
%
%   [W, T_out, X_dest, phase_out] = FCN_TURBINE( TURBINE, AMBIENT )
%
%   TURBINE fields: FLUID, mdot [kg/s], p_in, p_out [Pa], T_in [K], eta_is [-]
%   W          shaft power output [W]
%   T_out      outlet temperature [K]
%   X_dest     exergy destruction [W]
%   phase_out  CoolProp phase index of the exhaust (6 = two-phase, i.e. wet
%              expansion), evaluated from the actual (p, h) outlet state

% CoolProp failures (physically infeasible states) are caught so that the
% solver can recover; any other error is a genuine bug and is rethrown.
try
    [h_in, s_in] = fcn_props( TURBINE.FLUID, 'PT', TURBINE.p_in, TURBINE.T_in, 'H', 'S' );
    h_out_is = fcn_props( TURBINE.FLUID, 'PS', TURBINE.p_out, s_in, 'H' );
    h_out = h_in - TURBINE.eta_is*(h_in - h_out_is);
    if nargout > 3
        [T_out, s_out, phase_out] = fcn_props( TURBINE.FLUID, 'PH', TURBINE.p_out, h_out, 'T', 'S', 'PHASE' );
    else
        [T_out, s_out] = fcn_props( TURBINE.FLUID, 'PH', TURBINE.p_out, h_out, 'T', 'S' );
    end

    W = TURBINE.mdot*(h_in - h_out);
    X_dest = TURBINE.mdot*AMBIENT.T0*(s_out - s_in);
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    W = NaN;
    T_out = 0;
    X_dest = NaN;
    phase_out = NaN;
end

end
