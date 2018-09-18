function BuildBCMaps2D()

% function BuildMaps2DBC
% Purpose: Segregation of boundary face maps based on face tag

Globals2D_DG;

% create label of face nodes with boundary types from BCTag
ind                     = find(BCTag~=0);
bct                     = zeros(1,KG);
bct(1,EToGE(ind)')      = BCTag(ind)';
bnodes                  = ones(Nfp, 1)*bct(:)';
bnodes                  = bnodes(:);

% Create GEBC, map and vmap lists for BCs
GEBC_list   = containers.Map('KeyType','uint32','ValueType','any');
mapBC_list  = containers.Map('KeyType','uint32','ValueType','any');
vmapBC_list = containers.Map('KeyType','uint32','ValueType','any');

% find location of boundary nodes in face and volume node lists
for i=1:length(BC_ess_flags(:,1))
    tag    = BC_ess_flags(i,1);
    bc_ind = find(bnodes==tag);
    GEBC_list(tag)   = find(bct(1,:)==tag);
    mapBC_list(tag)  = reshape(mapB(bc_ind),Nfp,[]);
    vmapBC_list(tag) = reshape(vmapB(bc_ind),Nfp,[]);
end

return;
