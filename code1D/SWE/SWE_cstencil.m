function [q_lim] = SWE_cstencil(q,gravity,ind)

% Function converts conserved variable to characteristic variables
% locally on each 3-cell stencil, applies the limiter, and
% converts back to conserved variales

Globals1D_DG;

q_lim = q;

% Extend solution vector based on bc_type
depth_ext      = apply_bc(q(:,:,1));
discharge_ext  = apply_bc(q(:,:,2));

%depthc = zeros(1,K); dischargec = zeros(1,K); Char = zeros(2,Np,3);

% Compute cell averages
depthh = invV*depth_ext; depthh(2:Np,:)=0;
deptha = V*depthh; depthc = deptha(1,:);

dischageh = invV*discharge_ext; dischageh(2:Np,:)=0;
dischargea = V*dischageh; dischargec = dischargea(1,:);

% Compute characterisic variables
if(~isempty(ind))
    for i=1:length(ind)
        
        [L,invL] = SWECharMat(depthc(1,ind(i)+1),dischargec(1,ind(i)+1),gravity);
        
        Char  = invL*[depth_ext(:,ind(i))', depth_ext(:,ind(i)+1)', depth_ext(:,ind(i)+2)';
            discharge_ext(:,ind(i))', discharge_ext(:,ind(i)+1)', discharge_ext(:,ind(i)+2)'];
        
        Char1 = SlopeLimit3(reshape(Char(1,:,:),[Np,3]),x(:,ind(i)));
        Char2 = SlopeLimit3(reshape(Char(2,:,:),[Np,3]),x(:,ind(i)));
        
        U     = L*[Char1(:)';Char2(:)'];
        q_lim(:,ind(i),1) = reshape(U(1,:,:),[Np,1]);
        q_lim(:,ind(i),2) = reshape(U(2,:,:),[Np,1]);
        
    end
end


return;
