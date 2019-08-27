CleanUp1D;

clc
clear all
close all


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Sod';
rho_IC    =@(x) 1*(x<0.0) + 0.125*(x>=0.0);
vel_IC    =@(x) 0*x;
pre_IC    =@(x) 1*(x<0.0) + 0.1*(x>=0.0);

bnd_l     = -5;  
bnd_r     = 5;
mesh_pert = 0.0;
bc_cond   = {'D',1,'D',0.125;
             'D',0,'D',0.0;
             'D',1/(0.4),'D',0.1/(0.4)};  % For conserved variables
FinalTime = 2;
CFL       = 0.1;
K         = 200;
N         = 1;
RK        = 'LS54';


Indicator = 'NONE'; TVBM=1;
ind_var        = 'prim';
nn_model       = 'MLP_v1';	
Limiter    = 'NONE';
lim_var        = "char_stencil";

%Visc_model = 'NONE';
nn_visc_model = 'MLP_visc';
Visc_model='EV'; c_E=1; c_max=0.5;
%Visc_model='MDH'; c_A=2.5; c_k=0.2; c_max=0.5;
%Visc_model='MDA'; c_max=1;
%Visc_model='NN';
visc_var='density';

plot_iter  = 200;
save_iter  = 1;
save_soln  = true;
save_ind   = true;
save_visc  = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [0,1.2; -0.2,1.5; 0,1.2];

% Call code driver
EulerDriver1D; 




