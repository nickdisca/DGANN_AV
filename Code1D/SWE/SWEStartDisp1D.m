fprintf('Starting solve with DGANN for %s model, with the following parameters:\n',Problem.model)
fprintf('   test               : %s\n',Problem.test_name)
fprintf('   gravity            : %.2f\n',Problem.gravity)
fprintf('   mesh               : [%.3f,%.3f] with %d elements\n',Mesh.bnd_l,Mesh.bnd_r,Mesh.K)
fprintf('   FinalTime          : %.2f\n',Problem.FinalTime)
fprintf('   CFL                : %.2f\n',Problem.CFL)

if(strcmp(Indicator,'NONE') || strcmp(Indicator,'ALL') || strcmp(Indicator,'MINMOD'))
    fprintf('   Indicator          : %s \n',Limit.Indicator)
elseif(strcmp(Indicator,'NN'))
    fprintf('   Indicator          : %s (%s)\n',Limit.Indicator, Limit.nn_model)
elseif(strcmp(Indicator,'TVB'))
    fprintf('   Indicator          : %s (M=%.2f)\n',Limit.Indicator,Limit.TVBM)
end

if(~strcmp(Indicator,'NONE') && ~strcmp(Indicator,'ALL'))
    fprintf('   Ind_Var            : %s \n',Limit.ind_var)
end

fprintf('   Limiter            : %s\n',Limit.Limiter)

if(~strcmp(Limiter,'NONE'))
    fprintf('   Lim_Var            : %s \n',Limit.lim_var)
end