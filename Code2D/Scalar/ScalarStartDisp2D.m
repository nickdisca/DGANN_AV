fprintf('Starting solve with DGANN for %s model, with the following parameters:\n',Problem.model)
if(strcmp(Problem.model,'Advection'))
    fprintf('   test               : %s with velocity (%.2f,%.2f)\n',...
                                     Problem.test_name,Problem.AdvectionVelocity(1),Problem.AdvectionVelocity(2))
else
    fprintf('   test               : %s\n',Problem.test_name)
end   
fprintf('   mesh_file          : %s\n',Mesh.msh_file)
fprintf('   N                  : %d\n',Mesh.N)
if(~isempty(Problem.CFL))
    fprintf('   CFL                : %.2f\n',Problem.CFL)
else
    fprintf('   fixed_dt           : %.2e\n',Problem.fixed_dt)
end
fprintf('   tstamps            : %d\n',Problem.tstamps)

if(strcmp(Limit.Indicator,'NONE') || strcmp(Limit.Indicator,'ALL'))
    fprintf('   Indicator          : %s \n',Limit.Indicator)
elseif(strcmp(Limit.Indicator,'NN'))
    fprintf('   Indicator          : %s (%s)\n',Limit.Indicator, Limit.nn_model)
elseif(strcmp(Limit.Indicator,'TVB'))
    fprintf('   Indicator          : %s (M=%.2f, nu=%.2f)\n',Limit.Indicator,Limit.TVBM,Limit.TVBnu)
end

fprintf('   Limiter            : %s\n',Limit.Limiter)

switch Viscosity.model
    case 'NONE'
    fprintf('   Viscosity Model    : %s \n',Viscosity.model)
    
    case 'MDH'
    fprintf('   Viscosity Model    : %s (c_A=%.2f, c_k=%.2f, c_max=%.2f)\n',Viscosity.model,Viscosity.c_A,Viscosity.c_k,Viscosity.c_max)

    case 'MDA'
    fprintf('   Viscosity Model    : %s (c_max=%.2f)\n',Viscosity.model,Viscosity.c_max)
    
    case 'EV'
    fprintf('   Viscosity Model    : %s (c_E=%.2f, c_max=%.2f)\n',Viscosity.model,Viscosity.c_E,Viscosity.c_max)
    
    case 'NN'
    fprintf('   Viscosity Model    : %s (%s)\n',Viscosity.model, Viscosity.nn_visc_model)
end
