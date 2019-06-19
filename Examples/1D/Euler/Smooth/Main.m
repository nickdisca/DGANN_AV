CleanUp1D;

clc
clear all
close all


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Smooth';
rho_IC    =@(x) 1+0.5*sin(2*pi*x);
vel_IC    =@(x) ones(size(x));
pre_IC    =@(x) ones(size(x));

bnd_l     = 0;  
bnd_r     = 1;
mesh_pert = 0.0;
bc_cond   = {'P',0,'P',0;
             'P',0,'P',0.0;
             'P',0,'P',0};  % For conserved variables
FinalTime = 0.1;
CFL       = 0.1;
K         = 100;
N         = 1;
RK        = 'LS54';


Indicator = 'NONE'; TVBM=1;
ind_var        = 'prim';
nn_model       = 'MLP_v1';	
Limiter    = 'NONE';
lim_var        = "char_stencil";

nn_visc_model = 'MLP_visc';
%Visc_model = 'NONE';
%Visc_model='EV'; c_E=1; c_max=0.5;
%Visc_model='MDH'; c_A=2.5; c_k=0.2; c_max=0.5;
%Visc_model='MDA'; c_max=1;
Visc_model='NN';
visc_var='density';

plot_iter  = 20;
save_iter  = 1;
save_soln  = true;
save_ind   = true;
save_visc  = true;
save_plot  = true;
ref_avail  = false;
ref_fname  = 'ref_soln.dat';
var_ran    = [0,2; 0.5,1.5; 0.5,1.5];

% Call code driver
EulerDriver1D; 