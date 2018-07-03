CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Lax';
rho_IC =@(x) 0.445*(x<0) + 0.5*(x>=0.0);
vel_IC =@(x) 0.698*(x<0) + 0.0*(x>=0.0);
pre_IC =@(x) 3.528*(x<0) + 0.571*(x>=0.0);



bnd_l     = -5.0;  
bnd_r     = 5.0;
mesh_pert = 0.1;
bc_cond   = {'N',0,'N',0.0;
             'N',0,'N',0.0;
             'N',0,'N',0.0};
FinalTime = 1.3;
CFL       = 0.4;
Nelem     = 200;
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







