% Create BC_flag
CreateBC_Flags2D;

% Check parameters
ScalarCheckParam2D;


% Display paramaters
ScalarStartDisp2D;

% Initialize solver and construct grid and metric
[Mesh.VX,Mesh.VY,Mesh.K,Mesh.Nv,Mesh.EToV,Mesh.BFaces,Mesh.PerBToB_map,Mesh.PerBFToF_map] ...
                                          = read_gmsh_file(Mesh.msh_file);

% Generate necessary data structures
StartUp2D;

% Get essential BC_flags
Mesh.BC_ess_flags = BuildBCKeys2D(Mesh.BC_flags,Mesh.BC_ENUM.Periodic);

BuildBCMaps2D;

%% compute initial condition (time=0)
Q = feval(Problem.InitialCond, Mesh.x, Mesh.y);
    
% Find relative path
Find_relative_path;
    
% Extract MLP weights, biases and other parameters
if(strcmp(Limit.Indicator,'NN'))
    Net = read_mlp_param2D(Limit.nn_model,REL_PATH);
else
    Net.avail = false;
end
    
% Creating save file base names
Create_sfile_base2D;
    
% Solve Problem
fprintf('... starting main solve\n')

tic;
[Q_save,ind_save,ptc_hist,t_hist,Save_times] = Scalar2D(Q,Problem,Mesh,Limit,Net,Output);
sim_time = toc;
    
%%
% Saving data
if(Output.save_soln)
   Scalar_Save2D;
end
    
% Clean up processes
fprintf('... cleaning up\n')
CleanUp2D;

fprintf('------------ Solver has finished -------------\n')

