CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'LowDen';
rho_IC =@(x) 0*x + 7.0;
vel_IC =@(x) (x<0)*(-1.0)+ (x>=0.0)*1.0;
pre_IC =@(x) 0*x + 0.2;


bnd_l     = -1.0;  
bnd_r     = 1.0;
mesh_pert = 0.1;
bc_cond   = {'N',0,'N',0.0;
             'N',0,'N',0.0;
             'N',0,'N',0.0};
FinalTime = 0.6;
CFL       = 0.4;
Nelem     = 128;
N         = 4;


indicator_type = 'NN';
ind_var        = 'prim';
nn_model       = 'MLP_v1';
rec_limiter    = 'minmod';
lim_var        = "char_stencil";


plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
EulerDriver1D; 









