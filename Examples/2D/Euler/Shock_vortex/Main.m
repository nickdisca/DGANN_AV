%%

CleanUp2D;

close all
clear all
clc

%Model
model          = 'Euler';
gas_const      = 1.0;
gas_gamma      = 1.4;
test_name      = 'SV';
InitialCond    = @IC;
BC_cond        = {100001,'P'; 100002,'P'; 100003,'P'; 100004,'P'};

FinalTime      = 0.8;
CFL            = 0.4;
tstamps        = 1;
N              = 1;
RK               = 'LS54';

% Mesh file
msh_file      = 'nonuniform_m1p3_H002.msh';


% Set type of indicator
%Indicator = 'TVB'; TVBM = 200; TVBnu = 1.5;
Indicator     = 'NONE';
ind_var       = 'con';
nn_model      = 'MLP_v1';	
Limiter       = 'NONE'; 
lim_var       = 'con';
Filter_const  = true;


%Set viscosity model
%Visc_model = 'NONE';
Visc_model='EV'; c_E=1; c_max=0.25;
%Visc_model='MDH'; c_A=2; c_k=0.4; c_max=0.8;
%Visc_model='MDA'; c_max=0.8;
%Visc_model='NN';
visc_var = 'density';

% Mention which variables should be plotted
% Options available: 'density', 'velx', 'vely', 'pressure', 'energy'
plot_iter  = 50;
show_plot  = true;
xran       = [0,2]; 
yran       = [0,1];
plot_var   = {'density'};
clines     = {linspace(0.85,1.35,30)}; % Should have as many vectors as 
                                    % plotting variables
save_soln  = true;

    
EulerDriver2D;
