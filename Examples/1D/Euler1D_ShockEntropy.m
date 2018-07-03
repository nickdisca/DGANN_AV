CleanUp1D;

clc
clear all
close all


Globals1D_DG;
Globals1D_MLP;


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'ShockEntropy';
rho_IC =@(x) (x<-4)*3.857143 + (x>=-4).*(1 + 0.2*sin(5*x));
vel_IC =@(x) (x<-4)*2.629369;
pre_IC =@(x) (x<-4)*10.33333 + (x>=-4)*1.0;




bnd_l     = -5.0;  
bnd_r     = 5.0;
mesh_pert = 0.0;
bc_cond   = {'D',3.857143,'N',0.0;
             'D',10.141852,'D',0.0;
             'D',39.166661,'N',0.0};  % For conserved variables
FinalTime = 1.8;
CFL       = 0.1;
Nelem     = 200;
N         = 4;


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








