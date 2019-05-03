
% Check parameters
SWECheckParam1D;

% Display paramaters
SWEStartDisp1D;

% Find relative path
Find_relative_path;

% Generate simple mesh
[Mesh.Nv, Mesh.VX, Mesh.hK] = MeshGen1D(Mesh.bnd_l,Mesh.bnd_r,Mesh.K, Mesh.mesh_pert);

% generate various matrix operators and maps
StartUp1D;

% Extract MLP weights, biases and other parameters
if(strcmp(Limit.Indicator,'NN'))
    Net = read_mlp_param1D(Limit.nn_model,REL_PATH);
else
    Net.avail = false;
end

% Generate mass matrix and initialize solution
depth      = Problem.depth_IC(Mesh.x);
velocity   = Problem.vel_IC(Mesh.x);
discharge  = depth.*velocity;

% Creating vector of conserved variables
q = zeros(Mesh.Np,Mesh.K,2); q(:,:,1) = depth; q(:,:,2) = discharge;

% Creating file name base for saving solution
Output.fname_base = SWE_fnamebase1D(Problem,Mesh.N,Mesh.K,Limit,Mesh.mesh_pert);

if(Output.save_plot)
    assert(Output.save_soln && Output.save_ind,'To be able to save plots, set save_soln and save_ind as true');
end

% Solve Problem
fprintf('... starting main solve\n')
q   = SWE1D(q,Problem,Mesh,Limit,Net,Output);

% Save final solution
Save_SWE_soln1D(q,Mesh.x,Output.fname_base);

if(Output.save_plot)
    fprintf('... generating and saving plots in directory OUTPUT\n')
    if(Output.ref_avail)
        PlotSWE1D(Output.fname_base,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.rk_comb,true,Output.ref_fname); 
    else
        PlotSWE1D(Output.fname_base,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.rk_comb,false);
    end
end

% Clean up
fprintf('... cleaning up\n')
CleanUp1D;

fprintf('------------ Solver has finished -------------\n')
