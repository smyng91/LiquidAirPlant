function Ts = fcn_Ts( FLUID, p, T )
%FCN_TS State points for a T-s diagram.
%
%   Ts = FCN_TS( FLUID, p, T ) returns an N-by-2 matrix [s, T] for the N
%   states (p(i) [Pa], T(i) [K]); s is NaN where CoolProp cannot evaluate
%   the state. FLUID is a fluid string or a cell array with one per state.

T = T(:);
s = NaN(numel(T), 1);
for i = 1:numel(T)
    if iscell(FLUID), fluid = FLUID{i}; else, fluid = FLUID; end
    try
        s(i) = fcn_props( fluid, 'PT', p(i), T(i), 'S' );
    catch err
        if ~strcmp(err.identifier, 'fcn_props:failed'), rethrow(err); end
    end
end
Ts = [s, T];

end
