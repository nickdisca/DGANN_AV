%clc
%clear all
close all

% Global variables
Globals1D_DG;
Globals1D_MLP;

% Scalar Model
% model --> Advection, Burgers, BuckLev
model = 'BuckLev';

% Domain and time parameters by declaring variables
% Domain    --->  [bnd_l , bnd_r]
% bc_cond   ---> {'bc_type_left', bc_val_left, 'bc_type_right',bc_val_right}
%                  bc_types can be 'P', 'D' or 'N'
%                  bc_val is used only if bc_type is 'D'
% FinalTIme 
% CFL
% Nelem     ---> Number of cell/elements in the mesh
bnd_l     = 0.0;  
bnd_r     = 1.5;
mesh_pert = 0.1;
bc_cond   = {'N',0.0,'N',0.0};
FinalTime = 0.4;
CFL       = 0.4;
Nelem     = 150;

% Initial condition
test_name = 'Shockrare';
u_IC =@(x)  0.95*(x<0.5) +0.1*(x>=0.5);

% Order of polymomials used for approximation 
%N = 2;

% Troubled-cell indicator
% inidcator_type ---> minmod, TVB, NN
% TVB_M          ---> Parameter needed by TVB limiter
indicator_type = 'minmod';
indicator_type = 'TVB'; TVB_M = 10;
indicator_type = 'TVB'; TVB_M = 100;
indicator_type = 'TVB'; TVB_M = 1000;
indicator_type = 'NN';

% Neural Network Parameters
nn_model      = 'MLP_v1';

% Limiter used for reconstruction (this need not be the same as the 
% troubled-cell indicator)
% rec_limiter ---> none, minmod
rec_limiter = 'minmod';

% Plot and save parameters
plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
ScalarDriver1D; 




