CleanUp1D;

clc
clear all
close all


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
bc_cond   = {'D',1,'D',1.0;
             'D',0,'N',0.0;
             'D',1000/0.4,'D',0.01/0.4};
FinalTime = 0.012;
CFL       = 0.4;
K     = 256;
N         = 1;


Indicator = 'NN';
ind_var        = 'prim';
nn_model       = 'MLP_v1';
Limiter    = 'MINMOD';
lim_var        = "char_stencil";


% Plot and save parameters
plot_iter  = 10;
save_soln  = true;
save_ind   = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
rk_comb    = true;
var_ran    = [0,6; -1,21; 0,1100];

% Call code driver
EulerDriver1D; 





