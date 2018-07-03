CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


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
FinalTime = 0.4;
CFL       = 0.4;
Nelem     = 200;
N         = 3;


indicator_type = 'NN';
nn_model       = 'MLP_v1';
rec_limiter    = 'minmod';



plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
ScalarDriver1D; 




