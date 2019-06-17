% Check parameters
ScalarCheckParam1D;

% Display paramaters
ScalarStartDisp1D;

% Find relative path
Find_relative_path;

% Generate simple mesh
[Mesh.Nv, Mesh.VX, Mesh.hK] = MeshGen1D(Mesh.bnd_l,Mesh.bnd_r,Mesh.K, Mesh.mesh_pert);

% generate various matrix operators and maps
StartUp1D;

% Extract MLP weights, biases and other parameters for troubled cells
if(strcmp(Limit.Indicator,'NN'))
    Net = read_mlp_param1D(Limit.nn_model,REL_PATH);
else
    Net.avail = false;
end

%Repeat for viscosity
if(strcmp(Viscosity.model,'NN'))
    NetVisc = read_mlp_param1D_visc(Viscosity.nn_visc_model,REL_PATH,Mesh.N);
else
    NetVisc.avail = false;
end


% Generate mass matrix and initialize solution
u = Problem.u_IC(Mesh.x);

% Creating file name base for saving solution
Output.fname_base = Scalar_fnamebase1D(Problem,Mesh.N,Mesh.K,Limit,Viscosity,Mesh.mesh_pert);


if(save_plot)
    assert(save_soln && save_ind && save_visc,'To be able to save plots, set save_soln, save_ind, save_visc as true');
end

% Solve Problem
fprintf('... starting main solve\n')
[u] = Scalar1D(u,Problem,Mesh,Limit,Net,Viscosity,NetVisc,Output);

%%

% Save final solution
Save_scalar_soln1D(u,Mesh.x,Output.fname_base);

if(Output.save_plot)
    fprintf('... generating and saving plots in directory OUTPUT\n')
    if(Output.ref_avail)
        PlotScalar1D(Output.fname_base,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],true,Output.ref_fname,Mesh.x,Problem.RK); 
    else
        PlotScalar1D(Output.fname_base,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],false,Output.ref_fname,Mesh.x,Problem.RK);
    end
end

% Clean up
fprintf('... cleaning up\n')
CleanUp1D;

fprintf('------------ Solver has finished -------------\n')


