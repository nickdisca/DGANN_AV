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
CFL       = 0.6;
K     = 100;
N         = 1;


Indicator = 'MINMOD';
nn_model       = 'MLP_v1';	
Limiter    = 'MINMOD';


plot_iter  = 20;
save_soln  = true;
save_ind   = false;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
rk_comb    = false;
var_ran    = [-1.2,1.5];

% Call code driver
ScalarDriver1D; 




