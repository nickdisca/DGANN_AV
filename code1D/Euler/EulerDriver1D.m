% Find relative path
Find_relative_path;

% Generate simple mesh
[Nv, VX, K, hK] = MeshGen1D(bnd_l,bnd_r,Nelem, mesh_pert);

% generate various matrix operators and maps
StartUp1D;

% Extract MLP weights, biases and other parameters
if(strcmp(indicator_type,'NN'))
    [n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES,NN_Dir] = ...
        read_mlp_param1D(nn_model,REL_PATH);
end

% Assert BC condition is valid
Check_BC1D(bc_cond,3);

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
ind_fname = Euler_ind_fname1D(REL_PATH);

% Solve Problem
q = Euler1D(q,gas_gamma, gas_const,ind_fname);

% Save final solution
Save_Euler_soln1D(q,gas_gamma, gas_const,REL_PATH);

% Clean up
CleanUp1D;
