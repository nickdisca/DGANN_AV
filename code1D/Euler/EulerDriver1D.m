% Generate simple mesh
[Nv, VX, K, hK] = MeshGen1D(bnd_l,bnd_r,Nelem, mesh_pert);

% generate various matrix operators and maps
GenOps;

% Extract MLP weights, biases and other parameters
if(strcmp(indicator_type,'NN'))
    [n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES] = ...
        read_mlp_param(nn_model,sub_model,data_set,data_subset);
end

% Generate mass matrix and initialize solution
cx = ones(Np,1)*sum(M*x,1)/2;  
% depth = depth_IC(cx);
% vel   = vel_IC(cx);
rho        = rho_IC(x);
vel        = vel_IC(x);
pre        = pre_IC(x);
mmt        = rho.*vel;
energy     = 0.5*rho.*vel.^2 + pre/(gas_gamma - 1);

% Creating vector of conserved variables
q = zeros(Np,K,3); q(:,:,1) = rho; q(:,:,2) = mmt;  q(:,:,3) = energy;

% Creating file name for storing troubled-cell markers
ind_fname = Euler_ind_fname();

% Solve Problem
q = Euler1D(q,gas_gamma, gas_const,ind_fname);

% Save final solution
Save_Euler_soln(q,gas_gamma, gas_const);