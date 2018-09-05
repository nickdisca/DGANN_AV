function BuildBCMaps2D(BC_flags)

% function BuildMaps2DBC
% Purpose: Build specialized nodal maps for various types of
%          boundary conditions, specified in BCTag. 

Globals2D_DG;

% create label of face nodes with boundary types from BCTag
bct    = BCTag';
bnodes = ones(Nfp, 1)*bct(:)';
bnodes = bnodes(:);

% Create map and vmap lists for BCs
mapBC_list  = containers.Map('KeyType','uint32','ValueType','any');
vmapBC_list = containers.Map('KeyType','uint32','ValueType','any');

% find location of boundary nodes in face and volume node lists
for i=1:length(BC_flags(:,1))
    tag = BC_flags(i,1);
    if(BC_flags(i,2) ~= Periodic)
        mapBC = find(bnodes==tag);
        mapBC_list(tag)  = mapBC;
        vmapBC_list(tag) = vmapM(mapBC);
    end

end

return;
