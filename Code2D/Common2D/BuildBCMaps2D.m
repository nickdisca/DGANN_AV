function BuildBCMaps2D()

% function BuildMaps2DBC
% Purpose: Build specialized nodal maps for various types of
%          boundary conditions, specified in BCType. 

Globals2D_DG;

% create label of face nodes with boundary types from BCType
bct    = BCType';
bnodes = ones(Nfp, 1)*bct(:)';
bnodes = bnodes(:);

% find location of boundary nodes in face and volume node lists
mapI = find(bnodes==In);           vmapI = vmapM(mapI);
mapO = find(bnodes==Out);          vmapO = vmapM(mapO);
mapF = find(bnodes==Far);          vmapF = vmapM(mapF);
mapD = find(bnodes==Dirichlet);    vmapD = vmapM(mapD);
mapN = find(bnodes==Neuman);       vmapN = vmapM(mapN);
mapS = find(bnodes==Slip);         vmapS = vmapM(mapS);
return;
