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
CFL              = 0.3;  
tstamps          = 2;
N                = 1;
RK               = 'LS54';

% Set type of indicator
%Indicator       = 'TVB'; TVBM = 10; TVBnu = 1.5;
Indicator       = 'NONE';
Filter_const    = true;
nn_model        = 'MLP_v1';
Limiter         = 'NONE';


%Set viscosity model
%Visc_model = 'NONE';
nn_visc_model = 'MLP_visc';
Visc_model='EV'; c_E=1; c_max=0.25;
%Visc_model='MDH'; c_A=2; c_k=0.4; c_max=0.8;
%Visc_model='MDA'; c_max=0.8;
%Visc_model='NN';

% Mesh file
msh_file        = 'square_trans.msh';

% Output flags
plot_iter  = 50;
show_plot  = true;
xran       = [-0.5,0.5]; 
yran       = [-0.5,0.5]; 
clines     = linspace(-1.1,0.9,30);
save_soln  = true;

% Call main driver
ScalarDriver2D;

