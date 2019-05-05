CleanUp1D;

clc
clear all
close all



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
bc_cond   = {'D',1,'D',0.4;
             'D',1,'N',0.4;
             'D',3,'D',2.7};
FinalTime = 0.2;
CFL       = 0.4;
K         = 100;
N         = 4;


Indicator      = 'NN';
ind_var        = 'prim';
nn_model       = 'MLP_v1';
Limiter        = 'MINMOD';
lim_var        = "char_stencil";


% Plot and save parameters
plot_iter  = 100;
save_soln  = true;
save_ind   = true;
save_plot  = true;
ref_avail  = false;
ref_fname  = 'ref_soln.dat';
rk_comb    = true;
var_ran    = [0,1.2; 0.8,1.2; 0.8,1.2];

% Call code driver
EulerDriver1D; 






