CleanUp1D;

clc
clear all
close all

Globals1D_DG;
Globals1D_MLP;


model     = 'Advection';
test_name = 'Sine';
u_IC =@(x) sin(10*pi*x); 


bnd_l     = 0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 1.0;
CFL       = 0.2;
Nelem     = 100;
N         = 1;


indicator_type = 'minmod';
nn_model       = 'MLP_v1';	
rec_limiter    = 'minmod';


plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
ScalarDriver1D; 




