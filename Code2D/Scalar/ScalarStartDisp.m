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

fprintf('   Limiter            : %s\n',Limiter)
fprintf('   Mesf file          : %s\n\n',msh_file);