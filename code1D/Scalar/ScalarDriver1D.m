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
%u = u_IC(cx);
u = u_IC(x);

% Creating file name for storing troubled-cell markers
ind_fname = Scalar_ind_fname();

% Solve Problem
[u] = Scalar1D(u,ind_fname);

% Save final solution
Save_scalar_soln(u);




