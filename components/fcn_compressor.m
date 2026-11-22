function [W, T_out, X_dest, h_out] = fcn_compressor( COMPRESSOR, AMBIENT )
%FCN_COMPRESSOR Adiabatic compressor with isentropic efficiency.
%
%   [W, T_out, X_dest, h_out] = FCN_COMPRESSOR( COMPRESSOR, AMBIENT )
%
%   COMPRESSOR fields: FLUID, mdot [kg/s], p_in, p_out [Pa], T_in [K], eta_is [-]
%   W       shaft power input [W]
%   T_out   outlet temperature [K]
%   X_dest  exergy destruction [W]
%   h_out   outlet specific enthalpy [J/kg]

% CoolProp failures (physically infeasible states) are caught so that the
% solver can recover; any other error is a genuine bug and is rethrown.
try
    [h_in, s_in] = fcn_props( COMPRESSOR.FLUID, 'PT', COMPRESSOR.p_in, COMPRESSOR.T_in, 'H', 'S' );
    h_out_is = fcn_props( COMPRESSOR.FLUID, 'PS', COMPRESSOR.p_out, s_in, 'H' );
    h_out = h_in + (h_out_is - h_in)/COMPRESSOR.eta_is;
    [T_out, s_out] = fcn_props( COMPRESSOR.FLUID, 'PH', COMPRESSOR.p_out, h_out, 'T', 'S' );

    W = COMPRESSOR.mdot*(h_out - h_in);
    X_dest = COMPRESSOR.mdot*AMBIENT.T0*(s_out - s_in);
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    W = NaN;
    T_out = 0;
    X_dest = NaN;
    h_out = NaN;
end

end
