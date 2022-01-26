function [ f, g ] = fcn_evaluation( x, AMBIENT, PLANT, NTU, PTC, SHX )

PLANT.mdot_la = x(1);
PLANT.mdot_preheater = x(2);
PLANT.alpha_p_int = x(3);
PLANT.p_int = AMBIENT.p0 + (PLANT.p_pump-AMBIENT.p0)*PLANT.alpha_p_int;
PTC.mdot = x(4);
SHX.alpha1 = x(5);
SHX.mdot1 = PTC.mdot*SHX.alpha1;
SHX.alpha2 = 1-SHX.alpha1;
SHX.mdot2 = PTC.mdot*SHX.alpha2;
NTU.preheat = x(6);
NTU.recup = x(7);
NTU.shx1 = x(8);
NTU.shx2 = x(9);

% AAS config
options = optimoptions('fsolve','Display','None','UseParallel',false);
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );

if exitflag <= 0 || AAS.eta_II <= 0.001 || AAS.eta_rt <= 0.001
    AAS.eta_II = NaN;
    AAS.eta_rt = NaN;
end

% Objective functions F(X)
f(1) = -AAS.eta_rt;
f(2) = -AAS.eta_II;

% % Equality constraints G(X) = 0 MUST COME FIRST in g(1:me)
g = x(6)+x(7)+x(8)+x(9) - NTU.total;

% % Inequality constraints G(X) >= 0 MUST COME SECOND in g(me+1:m)
%   g(2) = x(2) - 1.333333333;
%   g(3) = x(3) - 2.666666666;

end