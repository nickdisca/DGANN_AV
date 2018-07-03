CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


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
Nelem     = 100;
N         = 1;


indicator_type = 'NN';
ind_var        = 'prim';
nn_model       = 'MLP_v1'; 
rec_limiter    = 'minmod';
lim_var        = 'char_stencil';



plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
SWEDriver1D; 




