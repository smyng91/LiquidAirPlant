clear; clc; close all

addpath('./components/');


%% optimization
% import design & operational parameters
plant_input;

%% setup MIDACO
key = 'Sam_Yang_(CAPS_Florida_State_Uni)_____[ACADEMIC-SINGLE-USER]';
problem.func = @(x) fcn_evaluation(x,AMBIENT, PLANT, NTU, PTC, SHX); % Call is [f,g] = problem_function(x)     
problem.o  = 2; % Number of objectives
problem.n  = 9; % Number of variables (in total)
problem.ni = 0; % Number of integer variables (0 <= nint <= n)
problem.m  = 1; % Number of constraints (in total)
problem.me = 1; % Number of equality constraints (0 <= me <=

problem.xl = [0.01, 0.01, 0, 0.1, 0.01, 0.1, 0.1, 0.1, 0.1 ];
problem.xu = [1, 10, 1, 20, 0.99, NTU.total, NTU.total, NTU.total, NTU.total ];
problem.x  = problem.xu; % Here for example: 'x' = lower bounds 'xl'
option.maxeval  = 100000;    % Maximum number of function evaluation (e.g. 1000000)
option.maxtime  = 60*60*24; % Maximum time limit in Seconds (e.g. 1 Day = 60*60*24)
option.printeval  = 1000;  % Print-Frequency for current best solution (e.g. 1000)
option.save2file  = 6000;     % Save SCREEN and SOLUTION to TXT-files [ 0=NO/ 1=YES]

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

