% Remove NN diretory paths if they still exist
CleanUp2D;

close all
clear all
clc

model             = 'Advection';
AdvectionVelocity = [1,1]; % Used for linear advection only
test_name         = 'Smooth'; 
InitialCond       = @IC;
BC_cond           = {100001,'P'; 100002,'P'; 100003,'P'; 100004,'P'};


FinalTime        = 0.5;
CFL              = 0.6;
tstamps          = 2;
N                = 1;

% Set type of indicator
Indicator       = 'TVB'; TVBM = 10; TVBnu = 1.5;
Indicator       = 'NN';
Filter_const    = true;
nn_model        = 'MLP_v1';
Limiter         = 'BJES';


% Mesh file
msh_file        = 'square_trans.msh';

% Output flags
plot_iter  = 50;
show_plot  = true;
xran       = [0,1]; 
yran       = [0,1];
clines     = linspace(-1,1,30);
save_soln  = true;

% Call main driver
ScalarDriver2D;

