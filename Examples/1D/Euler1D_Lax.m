%clc
%clear all
close all

% Global variables
Globals1D_DG;
Globals1D_MLP;

model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;

% Domain and time parameters by declaring variables
% Domain    --->  [bnd_l , bnd_r]
% bc_cond   ---> {'r_bc_type_left', r_bc_val_left, 'r_bc_type_right',r_bc_val_right;
%                 'm_bc_type_left', m_bc_val_left, 'm_bc_type_right',m_bc_val_right;
%                 'e_bc_type_left', e_bc_val_left, 'e_bc_type_right',e_bc_val_right}
%                  bc_types can be 'P', 'D' or 'N'
% FinalTIme 
% CFL
% Nelem     ---> Number of cell/elements in the mesh
bnd_l     = -5.0;  
bnd_r     = 5.0;
mesh_pert = 0.1;
bc_cond   = {'N',0,'N',0.0;
             'N',0,'N',0.0;
             'N',0,'N',0.0};
FinalTime = 1.3;
CFL       = 0.4;
Nelem     = 200;

% Initial condition
test_name = 'Lax';
rho_IC =@(x) 0.445*(x<0) + 0.5*(x>=0.0);
vel_IC =@(x) 0.698*(x<0) + 0.0*(x>=0.0);
pre_IC =@(x) 3.528*(x<0) + 0.571*(x>=0.0);

% Order of polymomials used for approximation 
N = 4;

% Troubled-cell indicator
% inidcator_type ---> minmod, TVB, NN
% TVB_M          ---> Parameter needed by TVB limiter
%indicator_type = 'minmod';
%indicator_type = 'TVB'; TVB_M = 10;
%indicator_type = 'TVB'; TVB_M = 100;
indicator_type = 'TVB'; TVB_M = 1000;
indicator_type = 'NN';

% Indicator variable
% ind_var ---> density
%              velocity
%              pressure
%              energy
%              prim
%              de
%ind_var = 'density';
%ind_var = 'velocity';
%ind_var = 'presure';
ind_var = 'prim';
%ind_var = 'energy';
%ind_var = 'de';


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
lim_var = "char_stencil";

% Plot and save parameters
plot_iter = 100;
save_soln = true;
save_ind  = true;

% Call code driver
EulerDriver1D; 




