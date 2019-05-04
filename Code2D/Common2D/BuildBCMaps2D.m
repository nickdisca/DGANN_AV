%function BuildBCMaps2D()

% function BuildMaps2DBC
% Purpose: Segregation of boundary face maps based on face tag

% create label of face nodes with boundary types from BCTag
ind                     = find(Mesh.BCTag~=0);
bct                     = zeros(1,Mesh.KG);
bct(1,Mesh.EToGE(ind)') = Mesh.BCTag(ind)';
bnodes                  = ones(Mesh.Nfp, 1)*bct(:)';
bnodes                  = bnodes(:);

% Create GEBC, map and vmap lists for BCs
Mesh.GEBC_list   = containers.Map('KeyType','uint32','ValueType','any');
Mesh.mapBC_list  = containers.Map('KeyType','uint32','ValueType','any');
Mesh.vmapBC_list = containers.Map('KeyType','uint32','ValueType','any');

% find location of boundary nodes in face and volume node lists
for i=1:length(Mesh.BC_ess_flags(:,1))
    tag    = Mesh.BC_ess_flags(i,1);
    bc_ind = find(bnodes==tag);
    Mesh.GEBC_list(tag)   = find(bct(1,:)==tag);
    Mesh.mapBC_list(tag)  = reshape(Mesh.mapB(bc_ind),Mesh.Nfp,[]);
    Mesh.vmapBC_list(tag) = reshape(Mesh.vmapB(bc_ind),Mesh.Nfp,[]);
end

%return;
