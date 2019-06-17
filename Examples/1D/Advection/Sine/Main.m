CleanUp1D;

clc
clear all
close all


model     = 'Advection';
test_name = 'Sine';
u_IC =@(x) sin(10*pi*x); 


bnd_l     = 0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 0.2;
CFL       = 0.2;
K     = 100;
N         = 3;
RK        = 'SSP3';


Indicator = 'NONE'; TVBM=1;
nn_model       = 'MLP_v1';	
Limiter    = 'NONE';

%Visc_model = 'NONE';
nn_visc_model = 'MLP_visc';
%Visc_model='EV'; c_E=1; c_max=0.5;
%Visc_model='MDH'; c_A=2.5; c_k=0.2; c_max=0.5;
%Visc_model='MDA'; c_max=1;
Visc_model='NN';

plot_iter  = 20;
save_soln  = true;
save_ind   = true;
save_visc  = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [-1.2,1.5];

% Call code driver
ScalarDriver1D; 




