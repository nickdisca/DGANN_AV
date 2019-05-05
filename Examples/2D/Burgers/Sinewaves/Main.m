% Remove NN diretory paths if they still exist
CleanUp2D;

close all
clear all
clc

model             = 'Burgers';
test_name         = 'Sinewaves'; 
InitialCond       = @IC;
BC_cond           = {100001,'P'; 100002,'P'; 100003,'P'; 100004,'P'};


FinalTime        = 0.3;
CFL              = 0.4;
%fixed_dt         = 1.0e-3;    
tstamps          = 2;
N                = 1;

% Set type of indicator
%Indicator       = 'TVB'; TVBM = 10; TVBnu = 1.5;
Indicator       = 'NN';
Filter_const    = true;
nn_model        = 'MLP_v1';
Limiter         = 'BJES';

% Mesh file
msh_file        = 'square.msh';

% Output flags
plot_iter  = 50;
show_plot  = true;
xran       = [-1,1]; 
yran       = [-1,1]; 
clines     = linspace(-1.2,1.2,30);
save_soln  = true;

% Call main driver
ScalarDriver2D;

