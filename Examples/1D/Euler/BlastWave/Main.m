CleanUp1D;

clc
clear all
close all


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Blast';
rho_IC =@(x) ones(size(x));
vel_IC =@(x) 0*x;
pre_IC =@(x) (x<0.1)*1000 + (x>=0.1).*(x<0.9)*0.01 + (x>=0.9)*100;


bnd_l     = 0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'N',0,'N',0.0;
             'D',0,'D',0.0;
             'N',0,'N',0.0};  
FinalTime = 0.038;
CFL       = 0.4;
K         = 400;
N         = 1;
RK        = 'LS54';


Indicator = 'NONE'; TVBM=1;
ind_var        = 'prim';
nn_model       = 'MLP_v1';	
Limiter    = 'NONE';
lim_var        = "char_stencil";

plot_iter  = 20;
save_iter  = 1;
save_soln  = true;
save_ind   = true;
save_visc  = true;
save_plot  = false;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [0,10; 0,4; 0,20];

nn_visc_model = 'MLP_visc';
%Visc_model = 'NONE';
visc_var='density';
nb_models=4;
Visc_model='EV'; c_E=1; c_max=0.5;


EulerDriver1D;