CleanUp1D;

clc
clear all
close all


model     = 'Burgers';
test_name = 'Complex';
u_IC =@(x)  sin(1*pi*x).*(abs(x)>=1)...
            + 3*(x>=-1).*(x<-0.5)...
            + 1*(x>=-0.5).*(x<0.0)...
            + 3*(x>=0.0).*(x<0.5)...
            + 2*(x>=0.5).*(x<1.0);


bnd_l     = -4.0;  
bnd_r     = 4.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 0.01;
CFL       = 0.4;
K     = 200;
N         = 2;
RK        = 'LS54';


Indicator = 'NONE'; TVBM=1;
nn_model       = 'MLP_v1';	
Limiter    = 'NONE';

nn_visc_model = 'MLP_visc';
Visc_model='EV'; c_E=1; c_max=0.5;


plot_iter  = 5;
save_iter  = 1;
save_soln  = true;
save_ind   = true;
save_visc  = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
rk_comb    = false;
var_ran    = [-1.2,4];

% Call code driver
ScalarDriver1D; 




