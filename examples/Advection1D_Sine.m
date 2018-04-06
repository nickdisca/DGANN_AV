clc
clear all
close all

% Global variables
Globals1D_DG;
Globals1D_MLP;

% Scalar Model
% model --> Advection, Burgers, BuckLev
model = 'Advection';

% Domain and time parameters by declaring variables
% Domain    --->  [bnd_l , bnd_r]
% mesh_pert --->  used to randomly perturbs uniform mesh by a factor of
%                 mesh_pert*(bnd_r-bnd_l)/Nelem. If this value is zero,
%                 a uniform mesh is generated
% bc_cond   ---> {'bc_type_left', bc_val_left, 'bc_type_right',bc_val_right}
%                  bc_types can be 'P', 'D' or 'N'
%                  bc_val is used only if bc_type is 'D'
% FinalTIme 
% CFL
% Nelem     ---> Number of cell/elements in the mesh
bnd_l     = 0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 1.0;
CFL       = 0.2;
Nelem     = 128;

% Initial condition
test_name = 'Sine';
% u_IC =@(x) sin(4*pi*x).*(x<0.5) + sin(8*pi*x).*(x<0.75).*(x>=0.5)...
%                              + sin(16*pi*x).*(x>=0.75);
u_IC =@(x) sin(10*pi*x);                         

% Order of polymomials used for approximation 
N = 4;

% Troubled-cell indicator
% inidcator_type ---> minmod, TVB, NN
% TVB_M          ---> Parameter needed by TVB limiter
indicator_type = 'minmod';
%indicator_type = 'TVB'; TVB_M = 10;
%indicator_type = 'TVB'; TVB_M = 100;
indicator_type = 'TVB'; TVB_M = 1000;
indicator_type = 'NN';

% Neural Network Parameters
nn_model      = 'MLP_v1';	

% Limiter used for reconstruction (this need not be the same as the 
% troubled-cell indicator)
% rec_limiter ---> none, minmod
% Note that if the TVB limiter is used as well, then the same TVB_M 
% parameter will be used for the reconstruction process
rec_limiter = 'minmod';

% Plot and save parameters
plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
ScalarDriver1D; 




