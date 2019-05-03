% Check parameters
EulerCheckParam1D;

% Display paramaters
EulerStartDisp1D;

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
rho        = Problem.rho_IC(Mesh.x);
vel        = Problem.vel_IC(Mesh.x);
pre        = Problem.pre_IC(Mesh.x);
mmt        = rho.*vel;
energy     = 0.5*rho.*vel.^2 + pre/(Problem.gas_gamma - 1);

% Creating vector of conserved variables
q = zeros(Mesh.Np,Mesh.K,3); q(:,:,1) = rho; q(:,:,2) = mmt;  q(:,:,3) = energy;

% Creating file name base for saving solution
Output.fname_base = Euler_fnamebase1D(Problem,Mesh.N,Mesh.K,Limit,Mesh.mesh_pert);


if(Output.save_plot)
    assert(Output.save_soln && Output.save_ind,'To be able to save plots, set save_soln and save_ind as true');
end

% Solve Problem
fprintf('... starting main solve\n')
q = Euler1D(q,Problem,Mesh,Limit,Net,Output);

% Save final solution
Save_Euler_soln1D(q,Mesh.x,Problem.gas_gamma,Output.fname_base);

%%
if(Output.save_plot)
    fprintf('... generating and saving plots in directory OUTPUT\n')
    if(Output.ref_avail)
        PlotEuler1D(Output.fname_base,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.rk_comb,true,Output.ref_fname); 
    else
        PlotEulerE1D(Output.fname_base,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.rk_comb,false);
    end
end

% Clean up
fprintf('... cleaning up\n')
CleanUp1D;

fprintf('------------ Solver has finished -------------\n')

