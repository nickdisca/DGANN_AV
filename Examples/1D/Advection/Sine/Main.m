CleanUp1D;

clc
clear all
close all


model     = 'Advection';
test_name = 'Sine';
u_IC =@(x) sin(10*pi*x); 


bnd_l     = 0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 2;
CFL       = 0.4;
K     = 100;
N         = 2;
RK        = 'SSP3';


Indicator = 'NONE';
nn_model       = 'MLP_v1';	
Limiter    = 'NONE';

Visc_model = 'NONE';
nn_visc_model = 'MLP_visc';


plot_iter  = 20;
save_soln  = true;
save_ind   = true;
save_visc  = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
rk_comb    = false;
var_ran    = [-1.2,1.5];

% Call code driver
ScalarDriver1D; 




