function ProjectFromNb2D = Get_Projection_2D
% Purspose: Build Projection maps (needed for Shu-Fu indicator)

Globals2D_DG;

ProjectFromNb2D  = zeros(Np,Np,3,K);

for i = 1:K
    for j = 1:3
        
        Enb = EToE(i,j);
        
        % Get vertex coordinates of neighbour jth neighbour
        v1 = EToV(Enb,1); x1 = VX(v1); y1 = VY(v1);
        v2 = EToV(Enb,2); x2 = VX(v2); y2 = VY(v2);
        v3 = EToV(Enb,3); x3 = VX(v3); y3 = VY(v3);
        
        % Correctional shift periodic ghost triangles
        if(isKey(PShift,i) && UseMeshPerData)
            pairdata = PShift(i);
            p_ind    = find(pairdata(:,1)==Enb);
            if ~isempty(p_ind)
                x1 = x1 - pairdata(p_ind,2);
                x2 = x2 - pairdata(p_ind,2);
                x3 = x3 - pairdata(p_ind,2);
                y1 = y1 - pairdata(p_ind,3);
                y2 = y2 - pairdata(p_ind,3);
                y3 = y3 - pairdata(p_ind,3);
            end
        end
        
        [re,se] = xy_ext_to_rs(x1,y1,x2,y2,x3,y3,x(:,i),y(:,i));
        
        
        
        % Find projection matrix
        ProjectFromNb2D(:,:,j,i)  = Vandermonde2D(N,re,se)*invV;
    end
end


return