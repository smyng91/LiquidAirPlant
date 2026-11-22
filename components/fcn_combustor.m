function [T_out, X_dest] = fcn_combustor( COMBUST, AMBIENT )
%FCN_COMBUSTOR Heat addition from natural gas combustion.
%
%   [T_out, X_dest] = FCN_COMBUSTOR( COMBUST, AMBIENT )
%
%   COMBUST fields: FLUID, mdot [kg/s], T_in [K], p_in, p_out [Pa],
%                   Q [W] (fuel heat input), NTU [-]
%   T_out   outlet temperature [K]
%   X_dest  exergy destruction [W]
%
%   A fraction 1 - exp(-NTU) of the fuel heat input reaches the working
%   fluid. The outlet state follows from the enthalpy balance
%   h_out = h_in + (1 - exp(-NTU))*Q/mdot, so the result does not depend on
%   a constant-cp approximation over the (large) temperature rise.

% CoolProp failures (physically infeasible states) are caught so that the
% solver can recover; any other error is a genuine bug and is rethrown.
try
    [h_in, s_in] = fcn_props( COMBUST.FLUID, 'PT', COMBUST.p_in, COMBUST.T_in, 'H', 'S' );
    h_out = h_in + (1 - exp(-COMBUST.NTU))*COMBUST.Q/COMBUST.mdot;
    [T_out, s_out] = fcn_props( COMBUST.FLUID, 'PH', COMBUST.p_out, h_out, 'T', 'S' );

    T_avg = (COMBUST.T_in + T_out)/2;
    X_dest = COMBUST.mdot*( h_in - h_out - AMBIENT.T0*(s_in - s_out) ) + ...
        COMBUST.Q*(1 - AMBIENT.T0/T_avg);
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    T_out = 0;
    X_dest = NaN;
end

end
