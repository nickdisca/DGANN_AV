function BC_ess_flags = BuildBCKeys2D(BC_flags,Pflag)

% function BuildBCKeys2D
% Purpose: Extracts BC Keys and types for non-periodic boundaries

BC_ess_flags = BC_flags;
ind = 0;

for i=1:length(BC_flags(:,1))
    if(BC_flags(i,2) ~= Pflag)
        ind = ind + 1;
        BC_ess_flags(ind,:) = BC_flags(i,:);
    end
end

BC_ess_flags = BC_ess_flags(1:ind,:);


return;
