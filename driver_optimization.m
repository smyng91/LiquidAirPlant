clear; clc; close all

addpath('./components/');


%% optimization
% import design & operational parameters
plant_input;

%% setup MIDACO
key = 'Sam_Yang_(CAPS_Florida_State_Uni)_____[ACADEMIC-SINGLE-USER]';
problem.func = @(x) problem_function(x,AMBIENT, PLANT, NTU, PTC, SHX); % Call is [f,g] = problem_function(x)     
problem.o  = 3; % Number of objectives
problem.n  = 8; % Number of variables (in total)
problem.ni = 0; % Number of integer variables (0 <= nint <= n)
problem.m  = 1; % Number of constraints (in total)
problem.me = 1; % Number of equality constraints (0 <= me <=

problem.xl = [0.01, 0.01, 0, 0.1, 0.01, 0.1, 0.1, 0.1 ];
problem.xu = [1, 10, 1, 20, 0.99, NTU.total, NTU.total, NTU.total ];
problem.x  = problem.xu; % Here for example: 'x' = lower bounds 'xl'
option.maxeval  = 80000;    % Maximum number of function evaluation (e.g. 1000000)
option.maxtime  = 60*60*24; % Maximum time limit in Seconds (e.g. 1 Day = 60*60*24)
option.printeval  = 1000;  % Print-Frequency for current best solution (e.g. 1000)
option.save2file  = 4000;     % Save SCREEN and SOLUTION to TXT-files [ 0=NO/ 1=YES]

option.param( 1) =  0;  % ACCURACY
option.param( 2) =  0;  % SEED
option.param( 3) =  0;  % FSTOP
option.param( 4) =  0;  % ALGOSTOP
option.param( 5) =  0;  % EVALSTOP
option.param( 6) =  0;  % FOCUS
option.param( 7) =  0;  % ANTS
option.param( 8) =  0;  % KERNEL
option.param( 9) =  0;  % ORACLE
option.param(10) =  0;  % PARETOMAX
option.param(11) =  0;  % EPSILON
option.param(12) =  0;  % BALANCE
option.param(13) =  0;  % CHARACTER

option.parallel = 20;  % Serial: 0 or 1, Parallel: 2,3,4,5,6,7,8...
[ solution ] = midaco( problem, option, key);


function [ f, g ] = problem_function( x, AMBIENT, PLANT, NTU, PTC, SHX )

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
NTU.shx2 = NTU.total-x(6)-x(7)-x(8);

% AAS config
options = optimoptions('fsolve','Display','None','UseParallel',false);
func = @(T) model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );
x0 = ones(12,1);
[T, fval, exitflag, output] = fsolve(func, x0, options);
[~, AAS, PUMP, PUMP_SHX, PREHEATER, RECUPERATOR, ...
    SHX, SHX1, SHX2, TURBINE_HP, TURBINE_LP, PTC] = ...
    model_AAS( T, AMBIENT, PLANT, NTU, PTC, SHX );

if exitflag <= 0 || AAS.eta_II <= 0 || AAS.eta_rt <= 0
    AAS.eta_II = NaN;
    AAS.eta_rt = NaN;
    AAS.COP = NaN;
end

% Objective functions F(X)
f(1) = -AAS.eta_rt;
f(2) = -AAS.eta_II;
f(3) = -AAS.COP;

% % Equality constraints G(X) = 0 MUST COME FIRST in g(1:me)
g(1) = x(6)+x(7)+x(8)+NTU.shx2 - NTU.total;
%   g(1) = x(1) - 1;
% % Inequality constraints G(X) >= 0 MUST COME SECOND in g(me+1:m)
%   g(2) = x(2) - 1.333333333;
%   g(3) = x(3) - 2.666666666;

end
