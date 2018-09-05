% Check parameters
ScalarCheckParam2D;

% Display paramaters
ScalarStartDisp;

% Initialize solver and construct grid and metric
[VX,VY,K,Nv,EToV,BFaces,PerBToB_map,PerBFToF_map] = read_gmsh_file(msh_file);

% Generate necessary data structures
StartUp2D;

% Turning non-perioidic BC faces to BC to Neumann
BCType = Neuman*(not(EToE - (1:K)'*ones(1,3)));

% Generates geometric data needed by Indicators
GetGeomIndData2D;

Q = zeros(Np, K, 1);

BuildBCMaps2D;

%% compute initial condition (time=0)

Ind_List = {'TVB2','TVB2', 'TVB2','NN_modal_patch_Pwise'};
Mlist = [10, 100, 200, 0];

% Ind_List = {'TVB'};
% Mlist = [10];
% nn_model_List = {'MLP_v2'};

for ilv = 1:3
     close all
%   
    Indicator = Ind_List{ilv}
    TVBM      = Mlist(ilv)
     ScalarCheckParam2D;
    Q = feval(InitialCond, x, y);
    
    % Find relative path
    Find_relative_path;
    
    % Extract MLP weights, biases and other parameters
    if(strcmp(Indicator,'NN') || ...
       strcmp(Indicator,'NN_Pwise') || ...
       strcmp(Indicator,'NN_modal_Pwise') || ...
       strcmp(Indicator,'NN_modal_patch_Pwise'))
        [n_input,n_output,n_hidden_layer,leaky_alpha,WEIGHTS,BIASES,NN_Dir] = ...
            read_mlp_param2D(nn_model,REL_PATH);
    end
    
    % Creating save file base names
    Create_sfile_base2D;
    
    % Solve Problem
    fprintf('... starting main solve\n')
    tic;
    [Q_save,ind_save,ptc_hist,t_hist,Save_times] = Scalar2D(Q,AdvectionVelocity,Save_soln);
    sim_time = toc;
    
    % Saving data
    Scalar_Save2D;
    
    % Clean up processes
    fprintf('... cleaning up\n')
    CleanUp2D;
end

