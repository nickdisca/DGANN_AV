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
Check_BC1D(bc_cond,1);

% Generate mass matrix and initialize solution
u = u_IC(x);

% Creating file name for storing troubled-cell markers
ind_fname = Scalar_ind_fname1D(REL_PATH);

% Solve Problem
[u] = Scalar1D(u,ind_fname);


% Save final solution
Save_scalar_soln1D(u,REL_PATH);

% Clean up
CleanUp1D;




