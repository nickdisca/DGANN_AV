CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Sine';
rho_IC =@(x) 1+0.5*sin(10*pi*x);
vel_IC =@(x) 0*x + 1;
pre_IC =@(x) 0*x + 1;



bnd_l     = -1.0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0,'P',0.0;
             'P',0,'P',0.0;
             'P',0,'P',0.0};  % For conserved variables
FinalTime = 1;
CFL       = 0.4;
Nelem     = 128;
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




