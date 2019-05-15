% Clean up various processes
if(exist('Net.NN_Dir','var'))
    if(~isempty(Net.NN_Dir))
        if(exist('Limit.Indicator','var'))
            if(strcmp(Limit.Indicator,'NN'))
                rmpath(Net.NN_Dir);
            end
        end
    end
end


