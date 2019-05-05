CleanUp1D;

clc
clear all
close all


model     = 'SWE';
gravity   = 1.0;
test_name = 'Dambreak';
depth_IC     =@(x) 3*(x<0.0) + 1*(x>=0.0);
velocity_IC  =@(x) 0*x;


bnd_l     = -3.0;  
bnd_r     = 3.0;
mesh_pert = 0.0;
bc_cond   = {'N',0.0,'N',0.0;
             'N',0.0,'N',0.0};  % For conserved variables
FinalTime = 1;
CFL       = 0.4;
K         = 100;
N         = 4;


Indicator = 'NN';
ind_var   = 'con';
nn_model  = 'MLP_v1'; 
Limiter   = 'MINMOD';
lim_var   = 'char_stencil';



plot_iter  = 10;
save_soln  = true;
save_ind   = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
rk_comb    = true;
var_ran    = [1,3.2; 0,1];

% Call code driver
SWEDriver1D; 




