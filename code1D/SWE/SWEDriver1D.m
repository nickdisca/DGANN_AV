% Generate simple mesh
[Nv, VX, K, hK] = MeshGen1D(bnd_l,bnd_r,Nelem, mesh_pert);

% generate various matrix operators and maps
GenOps1D;

% Extract MLP weights, biases and other parameters
if(strcmp(indicator_type,'NN'))
    [n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES] = ...
        read_mlp_param1D(nn_model);
end

% Assert BC condition is valid
Check_BC1D(bc_cond,2);

% Generate mass matrix and initialize solution
cx = ones(Np,1)*sum(M*x,1)/2;  
% depth = depth_IC(cx);
% vel   = vel_IC(cx);
depth      = depth_IC(x);
velocity   = velocity_IC(x);
discharge  = depth.*velocity;

% Creating vector of conserved variables
q = zeros(Np,K,2); q(:,:,1) = depth; q(:,:,2) = discharge;

% Creating file name for storing troubled-cell markers
ind_fname = SWE_ind_fname1D();

% Solve Problem
q = SWE1D(q,gravity,ind_fname);

% Save final solution
Save_SWE_soln1D(q);
