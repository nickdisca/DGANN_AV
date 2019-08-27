clear all; close all; clc;

addpath(genpath('.'));

nn_model='MLP_visc';

%polynomial order
N=4; 

%load the network
Net1D=read_mlp_param1D_visc(nn_model,'./',N);

%mesh consisting of 100 elements, domain is [0,1]
Nelem=100;
h=1/Nelem*ones(1,Nelem);

%construct the x coordinate
VX = (0:Nelem)/Nelem;
rv = JacobiGL(0,0,N);
x = ones(N+1,1)*VX(1:Nelem) + 0.5*(rv+1)*(VX(2:Nelem+1)-VX(1:Nelem));

%solution value
u=rand(N+1,Nelem);

%maximum local wave speed
local_wave_sp=ones(1,Nelem);

%boundary conditions
bc_cond   = {'P',0.0,'P',0.0};

%connectivity matrix
VtoE = zeros(2,Nelem);
for j=1:Nelem
    VtoE(1,j) = (j-1)*(N+1)+1;
    VtoE(2,j) = j*(N+1);
end

%predict the elementwise artificial viscosity coefficient
mu_piece = visc_MLP1D(u, Net1D, u, h ,local_wave_sp, VtoE, bc_cond);

%smoothen the viscosity
mu=Scalar1D_smooth_viscosity(mu_piece,x);
        