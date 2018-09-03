fprintf('Starting solve with DGANN for %s model, with the following parameters:\n',model)
fprintf('   N                  : %d\n',N)
if(strcmp(Indicator,'NN') || ...
   strcmp(Indicator,'NN_Pwise') || ...
   strcmp(Indicator,'NN_modal_Pwise') || ...
   strcmp(Indicator,'NN_modal_patch_Pwise'))
    fprintf('   Indicator          : %s (%s)\n',Indicator, nn_model)
elseif(strcmp(Indicator,'TVB'))
    fprintf('   Indicator          : %s (M=%.2f, nu=%.2f)\n',Indicator,TVBM,TVBnu)
end
fprintf('   Indicator variable : %s\n',ind_var)
fprintf('   Limiter            : %s\n',Limiter)
if(~strcmp(Indicator,'NONE'))
    fprintf('   Limiter variable   : %s\n',lim_var)
end
fprintf('   Mesf file          : %s\n\n',msh_file);