CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'IsoCon';
rho_IC =@(x) 1*(x<0) + 0.4*(x>=0.0);
vel_IC =@(x) 1*(x<0) + 1*(x>=0.0);
pre_IC =@(x) 1*(x<0) + 1*(x>=0.0);



bnd_l     = -1.0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'N',1,'N',0.4;
             'N',1,'N',1.0;
             'N',1,'N',1.0};
FinalTime = 0.5;
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






