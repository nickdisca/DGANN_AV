function [q_lim] = SWE_cstencil1D(q,gravity,ind,Limiter,Mesh)

% Function converts conserved variable to characteristic variables
% locally on each 3-cell stencil, applies the limiter, and
% converts back to conserved variales
% q includes ghost cell values

%Globals1D_DG;

q_lim = q(:,2:Mesh.K+1,:);

% Compute cell averages
% depthh = invV*q(:,:,1); depthh(2:Np,:)=0;
% deptha = V*depthh; depthc = deptha(1,:);
depthc = Mesh.AVG1D*q(:,:,1);

% dischageh = invV*q(:,:,2); dischageh(2:Np,:)=0;
% dischargea = V*dischageh; dischargec = dischargea(1,:);
dischargec = Mesh.AVG1D*q(:,:,2);

% Compute characterisic variables
if(~isempty(ind))
    for i=1:length(ind)
        
        [L,invL] = SWECharMat1D(depthc(1,ind(i)+1),dischargec(1,ind(i)+1),gravity);
        
        Char  = invL*[q(:,ind(i),1)', q(:,ind(i)+1,1)', q(:,ind(i)+2,1)';
                      q(:,ind(i),2)', q(:,ind(i)+1,2)', q(:,ind(i)+2,2)'];
        
        Char1 = SlopeLimit3(reshape(Char(1,:,:),[Mesh.Np,3]),Mesh.x(:,ind(i)),Limiter,Mesh);
        Char2 = SlopeLimit3(reshape(Char(2,:,:),[Mesh.Np,3]),Mesh.x(:,ind(i)),Limiter,Mesh);
        
        U     = L*[Char1(:)';Char2(:)'];
        q_lim(:,ind(i),1) = reshape(U(1,:,:),[Mesh.Np,1]);
        q_lim(:,ind(i),2) = reshape(U(2,:,:),[Mesh.Np,1]);
        
    end
end


return;
