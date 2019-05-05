CleanUp1D;

clc
clear all
close all


model     = 'Burgers';
test_name = 'Shocks';
u_IC =@(x)  10*(x<0.2) + 6*(x>=0.2).*(x<0.4)...
            + 0*(x>=0.4).*(x<0.6) -4*(x>=0.6);           

        
bnd_l     = 0.0;  
bnd_r     = 1.0;
mesh_pert = 0.1;
bc_cond   = {'N',0.0,'N',0.0};
FinalTime = 0.1;
CFL       = 0.2;
K     = 200;
N         = 1;


Indicator = 'MINMOD';
nn_model       = 'MLP_v1';
Limiter    = 'MINMOD';


plot_iter  = 100;
save_soln  = true;
save_ind   = true;
save_plot  = true;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
rk_comb    = false;
var_ran    = [-4.2,10.5];

% Call code driver
ScalarDriver1D; 




