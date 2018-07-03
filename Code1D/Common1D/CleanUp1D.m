% Clean up various processes
if(exist('NN_Dir','var'))
    if(~isempty(NN_Dir))
        if(exist('Indicator','var'))
            if(strcmp(Indicator,'NN'))
                rmpath(NN_Dir);
            end
        end
    end
end


