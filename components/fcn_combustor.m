function [T_out, X_dest, mdot_out] = fcn_combustor( COMBUST, AMBIENT )
%FCN_COMBUSTOR Heat addition from natural gas combustion.
%
%   [T_out, X_dest, mdot_out] = FCN_COMBUSTOR( COMBUST, AMBIENT )
%
%   COMBUST fields: FLUID, mdot [kg/s] (working fluid entering), T_in [K],
%                   p_in, p_out [Pa], Q [W] (fuel heat input, mdot_fuel*LHV),
%                   NTU [-] (heat release factor 1 - exp(-NTU); Inf for
%                   complete, adiabatic combustion),
%                   mdot_fuel [kg/s] (optional, default 0),
%                   T_fuel [K] (optional fuel supply temperature, default T0)
%   T_out     outlet temperature [K]
%   X_dest    exergy destruction [W]
%   mdot_out  outlet mass flow rate, mdot + mdot_fuel [kg/s]
%
%   The fuel joins the stream and the combustion products are treated as
%   the working fluid. The outlet state follows from the energy balance
%     mdot_out*h_out = mdot*h_in + mdot_fuel*h(T_fuel) + (1 - exp(-NTU))*Q,
%   so it does not rely on a constant-cp approximation over the large
%   temperature rise.

mdot_fuel = 0;
if isfield(COMBUST, 'mdot_fuel'), mdot_fuel = COMBUST.mdot_fuel; end
T_fuel = AMBIENT.T0;
if isfield(COMBUST, 'T_fuel'), T_fuel = COMBUST.T_fuel; end
mdot_out = COMBUST.mdot + mdot_fuel;
T0 = AMBIENT.T0;

% CoolProp failures (physically infeasible states) are caught so that the
% solver can recover; any other error is a genuine bug and is rethrown.
try
    [h_in, s_in] = fcn_props( COMBUST.FLUID, 'PT', COMBUST.p_in, COMBUST.T_in, 'H', 'S' );
    if mdot_fuel > 0
        [h_f, s_f] = fcn_props( COMBUST.FLUID, 'PT', COMBUST.p_in, T_fuel, 'H', 'S' );
    else
        h_f = 0; s_f = 0;
    end
    h_out = ( COMBUST.mdot*h_in + mdot_fuel*h_f + (1 - exp(-COMBUST.NTU))*COMBUST.Q )/mdot_out;
    [T_out, s_out] = fcn_props( COMBUST.FLUID, 'PH', COMBUST.p_out, h_out, 'T', 'S' );

    T_avg = (COMBUST.T_in + T_out)/2;
    X_dest = COMBUST.mdot*(h_in - T0*s_in) + mdot_fuel*(h_f - T0*s_f) ...
        - mdot_out*(h_out - T0*s_out) + COMBUST.Q*(1 - T0/T_avg);
catch err
    if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    T_out = 0;
    X_dest = NaN;
end

end
