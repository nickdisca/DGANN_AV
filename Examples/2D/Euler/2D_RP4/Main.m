% Remove NN diretory paths if they still exist
CleanUp2D;

close all
clear all
clc

%Model
model          = 'Euler';
gas_const      = 1.0;
gas_gamma      = 1.4;
test_name      = 'RP4';
InitialCond    = @IC;
BC_cond        = {100001,'P'; 100002,'P'; 100003,'P'; 100004,'P'};

FinalTime      = 0.25;
CFL            = 0.4;
%fixed_dt       = 1e-4;
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
msh_file      = 'square_trans.msh';

% Mention which variables should be plotted
% Options available: 'density', 'velx', 'vely', 'pressure', 'energy'
plot_iter  = 50;
show_plot  = true;
xran       = [-0.5,0.5]; 
yran       = [-0.5,0.5]; 
plot_var   = {'density','pressure'};
clines     = {linspace(0.255,1.9,30),linspace(0.36,2.15,30)}; % Should have as many vectors as 
                                                              % plotting variables
save_soln  = true;                                    

% Call main driver
EulerDriver2D;

