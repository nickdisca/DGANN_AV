fprintf('Starting solve with DGANN for %s model, with the following parameters:\n',Problem.model)
fprintf('   test               : %s\n',Problem.test_name)
fprintf('   gas_const          : %.2f\n',Problem.gas_const)
fprintf('   gas_gamma          : %.2f\n',Problem.gas_gamma)
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

if(~strcmp(Limit.Indicator,'NONE') && ~strcmp(Limit.Indicator,'ALL'))
    fprintf('   Ind_Var            : %s \n',Limit.ind_var)
end

fprintf('   Limiter            : %s\n',Limit.Limiter)

if(~strcmp(Limit.Limiter,'NONE'))
    fprintf('   Lim_Var            : %s \n',Limit.lim_var)
end
