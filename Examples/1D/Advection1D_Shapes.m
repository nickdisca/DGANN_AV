CleanUp1D;


clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Advection';
test_name = 'Shapes';
u_IC =@(x) 10*(x-0.2).*(x>=0.2).*(x<0.3)... 
            + 10*(0.4-x).*(x>=.3).*(x<0.4)...
            + 1*(x>=.6).*(x<0.8)...
            + 100*(x-1.0).*(1.2-x).*(x>=1.0).*(x<1.2); 

        
bnd_l     = 0.0;  
bnd_r     = 1.4;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 1.4;
CFL       = 0.2;
Nelem     = 100;
N         = 4;


indicator_type = 'TVB'; TVB_M = 10;
nn_model       = 'MLP_v2';	
rec_limiter    = 'minmod';



plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
ScalarDriver1D; 




