function [T, fval, exitflag, output] = solve_model( func, x0, options, name )
%SOLVE_MODEL Solve a plant model with fsolve and check the convergence.
%
%   [T, fval, exitflag, output] = SOLVE_MODEL( FUNC, X0, OPTIONS, NAME )
%
%   Wraps fsolve and issues a warning when no steady state was found, so
%   that results evaluated at a non-converged point are not mistaken for a
%   valid solution.

[T, fval, exitflag, output] = fsolve( func, x0, options );

res = max(abs(fval));
if exitflag <= 0 || ~(res < 1e-6)
    warning('LiquidAirPlant:notConverged', ...
        ['%s: fsolve found no steady state (exitflag %d, max |residual| = %.2e ', ...
         'relative to T0). The reported performance is NOT a valid solution; ', ...
         'check the inputs for feasibility.'], name, exitflag, res);
end

end
