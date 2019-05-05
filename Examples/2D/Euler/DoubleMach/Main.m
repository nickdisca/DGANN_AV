% Remove NN diretory paths if they still exist
CleanUp2D;

close all
clear all
clc

%Model
model          = 'Euler';
gas_const      = 1.0;
gas_gamma      = 1.4;
test_name      = 'DoubleMach';
InitialCond    = @IC;
BC_cond        = {101,'I'; 102,'O'; 103,'O'; 104,'S'; 105,'D'};

FinalTime      = 0.2;
CFL            = 0.4;
%fixed_dt       = 2e-5;%1e-5;
%fixed_dt       = 8e-6;%2e-5;
tstamps        = 1;
N              = 1;

% Set type of indicator
%Indicator = 'TVB'; TVBM = 200; TVBnu = 1.5;
Indicator     = 'NN';
ind_var       = 'con';
nn_model      = 'MLP_v1';	
Limiter       = 'BJES'; 
lim_var       = 'con';
Filter_const  = true;


% Mesh file
msh_file      = 'square.msh';

% Mention which variables should be plotted
% Options available: 'density', 'velx', 'vely', 'pressure', 'energy'
plot_iter  = 50;
show_plot  = true;
xran       = [0,4]; 
yran       = [0,1];
plot_var   = {'density'};
clines     = {linspace(1.5,21.5,30)}; 
save_soln  = true;

% Call main driver
EulerDriver2D;


