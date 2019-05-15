function ulimit = SlopeLimit1D(u,ind,Limiter,Mesh)

% Purpose: Apply slopelimiter (Pi^N) to u assuming u an N'th order polynomial
% u is must include ghost cell values

%Globals1D_DG;

ulimit = u(:,2:end-1);

if strcmp(Limiter,'NONE')
    return
else
     % Check to see if any elements require limiting
    if(~isempty(ind))
        
        ind_km1 = ind; 
        ind_k   = ind+1;
        ind_kp1 = ind+2;
    
        % Getting modal values
        uh = Mesh.invV*u; 
        
        % getting cell averages
        %uavg = uh; uavg(2:Np,:)=0; uavg = V*uavg; v = uavg(1,:);
        v = Mesh.AVG1D*u;
        vk = v(ind_k); vkm1 = v(ind_km1); vkp1 = v(ind_kp1);
        
        % create piecewise linear solution for limiting on specified elements
        uh(3:end,:)=0; ul = Mesh.V*uh(:,ind_k);
        
        % apply slope limiter to selected elements
        if strcmp(Limiter,'MINMOD')
            ulimit(:,ind) = SlopeLimitLin(ul,Mesh.x(:,ind),vkm1,vk,vkp1,Mesh);
        else
            error('Limiter %s not available!!',Limiter)
        end
    end

end
return
