% function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D()
% Purpose: Connectivity and boundary tables in the K # of Np elements


% number volume nodes consecutively
nodeids      = reshape(1:Mesh.K*Mesh.Np, Mesh.Np, Mesh.K);
gnodeids     = reshape(1:Mesh.KG*Mesh.Np, Mesh.Np, Mesh.KG);
Mesh.vmapM   = zeros(Mesh.Nfp, Mesh.Nfaces, Mesh.K); 
Mesh.vmapP   = zeros(Mesh.Nfp, Mesh.Nfaces, Mesh.K); 
Mesh.mapM    = (1:Mesh.K*Mesh.Nfp*Mesh.Nfaces)';     
Mesh.mapP    = reshape(Mesh.mapM, Mesh.Nfp, Mesh.Nfaces, Mesh.K);
Mesh.mapB    = reshape((1:Mesh.KG*Mesh.Nfp)', Mesh.Nfp, Mesh.KG);  
Mesh.vmapB   = zeros(Mesh.Nfp, Mesh.KG);
 
% find index of face nodes with respect to volume node ordering
for k1=1:Mesh.K
  for f1=1:Mesh.Nfaces
    Mesh.vmapM(:,f1,k1) = nodeids(Mesh.Fmask(:,f1), k1);
  end
end

one = ones(1, Mesh.Nfp);
for k1=1:Mesh.K
  for f1=1:Mesh.Nfaces
    % find neighbor
    k2 = Mesh.EToE(k1,f1); f2 = Mesh.EToF(k1,f1); % f2 may be negative for periodic boundaries
    
    % reference length of edge
    v1 = Mesh.EToV(k1,f1); v2 = Mesh.EToV(k1, 1+mod(f1,Mesh.Nfaces));
    refd = sqrt( (Mesh.VX(v1)-Mesh.VX(v2))^2 + (Mesh.VY(v1)-Mesh.VY(v2))^2 );
    
    % find perpendicalar shift of faces, which is non-zero for periodic
    % face pairs
    xs = 0; ys = 0;
    if(isKey(Mesh.PShift,k1) && Mesh.UseMeshPerData)
        pairdata = Mesh.PShift(k1);
        p_ind    = find(pairdata(:,1)==k2);
        if ~isempty(p_ind)
            xs = pairdata(p_ind,2);
            ys = pairdata(p_ind,3);
        end
    end
    
    % find find volume node numbers of left and right nodes 
    vidM = Mesh.vmapM(:,f1,k1); vidP = Mesh.vmapM(:,f2,k2);
    x1 = Mesh.x(vidM); y1 = Mesh.y(vidM); x2 = Mesh.x(vidP)-xs; y2 = Mesh.y(vidP)-ys;
    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;
   
    % Compute distance matrix
    D = (x1 -x2').^2 + (y1-y2').^2;
    [idM, idP] = find(sqrt(abs(D))<Mesh.NODETOL*refd);
    assert(length(idM) == Mesh.Nfp);
    Mesh.vmapP(idM,f1,k1) = vidP(idP); 
    Mesh.mapP(idM,f1,k1) = idP + (f2-1)*Mesh.Nfp+(k2-1)*Mesh.Nfaces*Mesh.Nfp;
    
    % Find maps if triangle has a non-periodic boundary face
    % Note that the shared face is the third face of ghost element
    % vmapB stores loaction of boundary triangle face nodes that map
    % to the given ghost triangles third face. In otherwords, its a mapping
    % from ghost triangle faces to boundary triangle faces and NOT the
    % other way around.
    if(k2 == k1)
        kg2  = Mesh.EToGE(k1,f1);
        vidP = gnodeids(Mesh.Fmask(:,Mesh.Nfaces), kg2);
        x2   = Mesh.xG(vidP)-xs; y2 = Mesh.yG(vidP)-ys;
        x2   = x2*one;      y2 = y2*one;
        D = (x1 -x2').^2 + (y1-y2').^2;
        [idM, idP] = find(sqrt(abs(D))<Mesh.NODETOL*refd);
        assert(length(idM) == Mesh.Nfp);
        Mesh.vmapB(idP,kg2) = vidM(idM); 
        Mesh.mapB(idP,kg2) = idM + (f1-1)*Mesh.Nfp+(k1-1)*Mesh.Nfaces*Mesh.Nfp; 
    end
    
  end
end

% reshape vmapM and vmapP to be vectors and create boundary node list
Mesh.vmapP = Mesh.vmapP(:); Mesh.vmapM = Mesh.vmapM(:); Mesh.mapP = Mesh.mapP(:); 
Mesh.vmapM = reshape(Mesh.vmapM, Mesh.Nfp*Mesh.Nfaces, Mesh.K); 
Mesh.vmapP = reshape(Mesh.vmapP, Mesh.Nfp*Mesh.Nfaces, Mesh.K);

% nnx = nx(mapB); nnx = nnx(:);
% nny = ny(mapB); nny = nny(:);
% nnx = nnx./(sqrt(nnx.^2 + nny.^2));
% nny = nny./(sqrt(nnx.^2 + nny.^2));
% 
% figure(1000)
% plot(nnx,nny,'o')
% hold all



