% Clean up various processes
if(exist('NN_Dir','var'))
    if(~isempty(NN_Dir))
        if(exist('Indicator','var'))
            if(strcmp(Indicator,'NN') || ...
                    strcmp(Indicator,'NN_Pwise') || ...
                    strcmp(Indicator,'NN_modal_Pwise') || ...
                    strcmp(Indicator,'NN_modal_patch_Pwise'))
                rmpath(NN_Dir);
            end
        end
    end
end


