% Remove NN diretory paths if they still exist
CleanUp2D;

close all
clear all
clc

model             = 'Burgers';
test_name         = 'RP'; 
InitialCond       = @IC;
BC_cond           = {100001,'P'; 100002,'P'; 100003,'P'; 100004,'P'};


FinalTime        = 0.3;
fixed_dt         = 1.0e-3;    
tstamps          = 2;
N                = 3;

% Set type of indicator
%Indicator       = 'TVB'; TVBM = 10; TVBnu = 1.5;
Indicator       = 'NN';
Filter_const    = true;
nn_model        = 'MLP_v1';
Limiter         = 'BJES';


% Mesh file
msh_file        = 'square_trans.msh';

% Output flags
plot_iter  = 50;
show_plot  = true;
xran       = [-0.5,0.5]; 
yran       = [-0.5,0.5]; 
clines     = linspace(-0.98,0.78,30);
save_soln  = true;

% Call main driver
ScalarDriver2D;

