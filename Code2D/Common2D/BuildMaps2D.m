function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D()

% function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D
% Purpose: Connectivity and boundary tables in the K # of Np elements

Globals2D_DG;

% number volume nodes consecutively
nodeids = reshape(1:K*Np, Np, K);
gnodeids= reshape(1:KG*Np, Np, KG);
vmapM   = zeros(Nfp, Nfaces, K); vmapP   = zeros(Nfp, Nfaces, K); 
mapM    = (1:K*Nfp*Nfaces)';     mapP    = reshape(mapM, Nfp, Nfaces, K);
mapB   = reshape((1:KG*Nfp)', Nfp, KG);  vmapB   = zeros(Nfp, KG);
 
% find index of face nodes with respect to volume node ordering
for k1=1:K
  for f1=1:Nfaces
    vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
  end
end

one = ones(1, Nfp);
for k1=1:K
  for f1=1:Nfaces
    % find neighbor
    k2 = EToE(k1,f1); f2 = EToF(k1,f1); % f2 may be negative for periodic boundaries
    
    % reference length of edge
    v1 = EToV(k1,f1); v2 = EToV(k1, 1+mod(f1,Nfaces));
    refd = sqrt( (VX(v1)-VX(v2))^2 + (VY(v1)-VY(v2))^2 );
    
    % find perpendicalar shift of faces, which is non-zero for periodic
    % face pairs
    xs = 0; ys = 0;
    if(isKey(PShift,k1) && UseMeshPerData)
        pairdata = PShift(k1);
        p_ind    = find(pairdata(:,1)==k2);
        if ~isempty(p_ind)
            xs = pairdata(p_ind,2);
            ys = pairdata(p_ind,3);
        end
    end
    
    % find find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);
    x1 = x(vidM); y1 = y(vidM); x2 = x(vidP)-xs; y2 = y(vidP)-ys;
    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;
   
    % Compute distance matrix
    D = (x1 -x2').^2 + (y1-y2').^2;
    [idM, idP] = find(sqrt(abs(D))<NODETOL*refd);
    assert(length(idM) == Nfp);
    vmapP(idM,f1,k1) = vidP(idP); mapP(idM,f1,k1) = idP + (f2-1)*Nfp+(k2-1)*Nfaces*Nfp;
    
    % Find maps if triangle has a non-periodic boundary face
    % Note that the shared face is the third face of ghost element
    % vmapB stores loaction of boundary triangle face nodes that map
    % to the given ghost triangles third face. In otherwords, its a mapping
    % from ghost triangle faces to boundary triangle faces and NOT the
    % other way around.
    if(k2 == k1)
        kg2  = EToGE(k1,f1);
        vidP = gnodeids(Fmask(:,Nfaces), kg2);
        x2   = xG(vidP)-xs; y2 = yG(vidP)-ys;
        x2   = x2*one;      y2 = y2*one;
        D = (x1 -x2').^2 + (y1-y2').^2;
        [idM, idP] = find(sqrt(abs(D))<NODETOL*refd);
        assert(length(idM) == Nfp);
        vmapB(idP,kg2) = vidM(idM); mapB(idP,kg2) = idM + (f1-1)*Nfp+(k1-1)*Nfaces*Nfp; 
    end
    
  end
end

% reshape vmapM and vmapP to be vectors and create boundary node list
vmapP = vmapP(:); vmapM = vmapM(:); mapP = mapP(:); 
vmapM = reshape(vmapM, Nfp*Nfaces, K); vmapP = reshape(vmapP, Nfp*Nfaces, K);

% nnx = nx(mapB); nnx = nnx(:);
% nny = ny(mapB); nny = nny(:);
% nnx = nnx./(sqrt(nnx.^2 + nny.^2));
% nny = nny./(sqrt(nnx.^2 + nny.^2));
% 
% figure(1000)
% plot(nnx,nny,'o')
% hold all




return
