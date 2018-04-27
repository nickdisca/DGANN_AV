% clc
% clear all
close all

% Global variables
Globals1D_DG;
Globals1D_MLP;

model   = 'SWE';
gravity = 1.0;

% Domain and time parameters by declaring variables
% Domain    --->  [bnd_l , bnd_r]
% bc_cond   ---> {'dep_bc_type_left', dep_bc_val_left, 'dep_bc_type_right',dep_bc_val_right;
%                 'dis_bc_type_left', dis_bc_val_left, 'dis_bc_type_right',dis_bc_val_right }
%                  bc_types can be 'P', 'D' or 'N'
% FinalTIme 
% CFL
% Nelem     ---> Number of cell/elements in the mesh
bnd_l     = -3.0;  
bnd_r     = 3.0;
mesh_pert = 0.0;
bc_cond   = {'N',0.0,'N',0.0;
             'N',0.0,'N',0.0};  % For conserved variables
FinalTime = 1;
CFL       = 0.4;
Nelem     = 100;

% Initial condition
test_name = 'Dambreak';
depth_IC     =@(x) 3*(x<0.0) + 1*(x>=0.0);
velocity_IC  =@(x) 0*x;

% Order of polymomials used for approximation 
N = 4;

% Troubled-cell indicator
% inidcator_type ---> minmod, TVB, NN
% TVB_M          ---> Parameter needed by TVB limiter
indicator_type = 'minmod';
indicator_type = 'TVB'; TVB_M = 10;
%indicator_type = 'TVB'; TVB_M = 100;
%indicator_type = 'TVB'; TVB_M = 1000;
indicator_type = 'NN';
%indicator_type = 'FuShu';

% Indicator variable
% ind_var ---> depth
%              velocity
%              prim (depth and velocity) earlier set as dv
ind_var = 'depth';
ind_var = 'velocity';
ind_var = 'prim';


% Neural Network Parameters
nn_model      = 'MLP_v1'; 

% Limiter used for reconstruction (this need not be the same as the 
% troubled-cell indicator)
% rec_limiter ---> none, minmod
rec_limiter = 'minmod';

% Limiting variable
% lim_var ---> prim (primitive)
%              con (conserved)
%              char_cell (cell-wise transformed characteristic variables)
%              char_stencil (stencil-wise transformed characteristic variables)
%lim_var = "prim";
%lim_var = "con";
%lim_var = "char_cell";
lim_var = 'char_stencil';

% Plot and save parameters
plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
SWEDriver1D; 




