CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


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
Nelem     = 100;
N         = 1;


indicator_type = 'NN';
nn_model       = 'MLP_v1';
rec_limiter    = 'minmod';


plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
ScalarDriver1D; 




