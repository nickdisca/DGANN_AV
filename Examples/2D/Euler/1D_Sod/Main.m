% Remove NN diretory paths if they still exist
CleanUp2D;

close all
clear all
clc

%Model
model          = 'Euler';
gas_const      = 1.0;
gas_gamma      = 1.4;
test_name      = '1DSod';
InitialCond    = @IC;
BC_cond        = {101,'D'; 102,'D'; 103,'P'; 104,'P'};

FinalTime      = 1;
CFL            = 0.5;
%fixed_dt       = 4e-5;
N              = 3;
tstamps        = 4;


% Set type of indicator
Indicator = 'TVB'; TVBM = 10; TVBnu = 1.5;
Indicator     = 'NN';
ind_var       = 'con';
nn_model      = 'MLP_v1';	
Limiter       = 'BJES'; 
lim_var       = 'con';
Filter_const  = true;
Remove_iso    = false;

% Mesh file
msh_file      = 'square.msh';

% Mention which variables should be plotted
% Options available: 'density', 'velx', 'vely', 'pressure', 'energy'
plot_iter  = 50;
show_plot  = false;
xran       = [-3,3]; 
yran       = [-.1,.1]; 
plot_var   = {'density'};
clines     = {linspace(0.126,0.99,30)}; % Should have as many vectors as 
                                      % plotting variables
save_soln  = true;

% Call main driver
EulerDriver2D;

