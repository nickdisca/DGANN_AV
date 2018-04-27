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
Check_BC1D(bc_cond,2);

% Generate mass matrix and initialize solution
depth      = depth_IC(x);
velocity   = velocity_IC(x);
discharge  = depth.*velocity;

% Creating vector of conserved variables
q = zeros(Np,K,2); q(:,:,1) = depth; q(:,:,2) = discharge;

% Creating file name for storing troubled-cell markers
ind_fname = SWE_ind_fname1D(REL_PATH);

% Solve Problem
q = SWE1D(q,gravity,ind_fname);

% Save final solution
Save_SWE_soln1D(q,REL_PATH);

% Clean up
CleanUp1D;
