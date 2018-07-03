CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'LeftBlastWave';
rho_IC =@(x) 0*x + 1.0;
vel_IC =@(x) 0*x;
pre_IC =@(x) 1000.0*(x<0.5) + 0.01*(x>=0.5);




bnd_l     = 0;  
bnd_r     = 1;
mesh_pert = 0.1;
bc_cond   = {'N',0,'N',0.0;
             'N',0,'N',0.0;
             'N',0,'N',0.0};
FinalTime = 0.012;
CFL       = 0.4;
Nelem     = 256;
N         = 1;


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





