function varargout = fcn_props( fluid, inputs, v1, v2, varargin )
%FCN_PROPS Thermophysical properties from CoolProp, one flash per call.
%
%   [out1, out2, ...] = FCN_PROPS( FLUID, INPUTS, V1, V2, NAME1, NAME2, ... )
%
%   FLUID   CoolProp fluid string, e.g. 'AIR' or 'INCOMP::S800'
%   INPUTS  'PT'  V1 = pressure [Pa], V2 = temperature [K]
%           'PH'  V1 = pressure [Pa], V2 = specific enthalpy [J/kg]
%           'PS'  V1 = pressure [Pa], V2 = specific entropy [J/kg/K]
%   NAMEk   'T' [K], 'P' [Pa], 'H' [J/kg], 'S' [J/kg/K], 'CP' [J/kg/K],
%           'D' [kg/m3], 'MU' [Pa s], 'K' [W/m/K], 'PR' [-], 'PHASE' [-]
%           (PHASE uses CoolProp's phase indices, e.g. 6 = two-phase)
%
%   Example: [h, s] = fcn_props('AIR', 'PT', 1e5, 298, 'H', 'S');
%
%   A CoolProp AbstractState is created once per fluid and reused, so all
%   requested properties come from a single flash calculation. This is much
%   faster than calling PropsSI once per property. If the low-level
%   interface is unavailable, PropsSI is used instead (same results).
%
%   Any CoolProp failure (state outside the valid range, non-finite input,
%   etc.) is raised as an error with identifier 'fcn_props:failed', which the
%   component models catch to keep the solver alive in infeasible regions.

persistent states pair use_low_level

if isempty(use_low_level)
    try
        ref = double(py.CoolProp.CoolProp.PropsSI('Hmass','P',101325,'T',300,'AIR'));
    catch err
        error('fcn_props:noCoolProp', ...
            'CoolProp for Python is not reachable from MATLAB: %s', err.message);
    end
    states = containers.Map();
    try
        % input-pair codes differ between CoolProp versions: query them
        pair.PT = py.CoolProp.CoolProp.PT_INPUTS;
        pair.PH = py.CoolProp.CoolProp.HmassP_INPUTS;
        pair.PS = py.CoolProp.CoolProp.PSmass_INPUTS;
        % self-test: the low-level interface must reproduce PropsSI
        AS = new_state('AIR');
        AS.update(pair.PT, 101325, 300);
        use_low_level = abs(double(AS.hmass()) - ref) <= 1e-9*abs(ref);
        if use_low_level
            states('AIR') = AS;
        end
    catch
        use_low_level = false;
    end
    if ~use_low_level
        warning('fcn_props:fallback', ...
            'CoolProp AbstractState unavailable; falling back to (slower) PropsSI.');
    end
end

if numel(varargin) < max(nargout, 1)
    error('fcn_props:badRequest', 'Requested %d outputs but named %d properties.', ...
        max(nargout, 1), numel(varargin));
end
if ~(isreal(v1) && isreal(v2) && isfinite(v1) && isfinite(v2))
    error('fcn_props:failed', 'Non-finite or complex input for %s.', fluid);
end

varargout = cell(1, max(nargout, 1));
try
    if use_low_level
        if isKey(states, fluid)
            AS = states(fluid);
        else
            try
                AS = new_state(fluid);
            catch err
                error('fcn_props:badRequest', 'Unknown CoolProp fluid ''%s'': %s', ...
                    fluid, err.message);
            end
            states(fluid) = AS;
        end
        switch inputs
            case 'PT', AS.update(pair.PT, v1, v2);
            case 'PH', AS.update(pair.PH, v2, v1);   % CoolProp order: (h, p)
            case 'PS', AS.update(pair.PS, v1, v2);
            otherwise, bad_request(inputs);
        end
        for k = 1:numel(varargout)
            switch varargin{k}
                case 'T',     val = AS.T();
                case 'P',     val = AS.p();
                case 'H',     val = AS.hmass();
                case 'S',     val = AS.smass();
                case 'CP',    val = AS.cpmass();
                case 'D',     val = AS.rhomass();
                case 'MU',    val = AS.viscosity();
                case 'K',     val = AS.conductivity();
                case 'PR',    val = AS.Prandtl();
                case 'PHASE', val = py.int(AS.phase());
                otherwise,    bad_request(varargin{k});
            end
            varargout{k} = double(val);
        end
    else
        switch inputs
            case 'PT', n1 = 'P'; n2 = 'T';
            case 'PH', n1 = 'P'; n2 = 'Hmass';
            case 'PS', n1 = 'P'; n2 = 'Smass';
            otherwise, bad_request(inputs);
        end
        for k = 1:numel(varargout)
            varargout{k} = double(py.CoolProp.CoolProp.PropsSI( ...
                propssi_name(varargin{k}), n1, v1, n2, v2, fluid));
        end
    end
catch err
    if strcmp(err.identifier, 'fcn_props:badRequest')
        rethrow(err);
    end
    error('fcn_props:failed', '%s', err.message);
end

for k = 1:numel(varargout)
    if ~isfinite(varargout{k})
        error('fcn_props:failed', 'CoolProp returned a non-finite %s for %s.', ...
            varargin{k}, fluid);
    end
end

end

function AS = new_state( fluid )
% 'INCOMP::S800' -> backend 'INCOMP', fluid 'S800'; plain names use HEOS
idx = strfind(fluid, '::');
if isempty(idx)
    AS = py.CoolProp.CoolProp.AbstractState('HEOS', fluid);
else
    AS = py.CoolProp.CoolProp.AbstractState(fluid(1:idx(1)-1), fluid(idx(1)+2:end));
end
end

function name = propssi_name( key )
switch key
    case 'T',     name = 'T';
    case 'P',     name = 'P';
    case 'H',     name = 'Hmass';
    case 'S',     name = 'Smass';
    case 'CP',    name = 'CPMASS';
    case 'D',     name = 'Dmass';
    case 'MU',    name = 'VISCOSITY';
    case 'K',     name = 'CONDUCTIVITY';
    case 'PR',    name = 'PRANDTL';
    case 'PHASE', name = 'Phase';
    otherwise,    bad_request(key);
end
end

function bad_request( key )
error('fcn_props:badRequest', 'Unknown property or input pair ''%s''.', key);
end
