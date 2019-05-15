function [q_lim] = Euler_cstencil1D(q,gas_gamma,gas_const,ind,Limiter,Mesh)

% Function converts conserved variable to characteristic variables
% locally on each 3-cell stencil, applies the limiter, and
% converts back to conserved variales
% q includes ghost cell values

q_lim = q(:,2:Mesh.K+1,:);

% Extend solution vector based on bc_type
rho_ext   = q(:,:,1);
mmt_ext   = q(:,:,2);
Ener_ext  = q(:,:,3);

%depthc = zeros(1,K); dischargec = zeros(1,K); Char = zeros(2,Np,3);

% Compute cell averages
% rhoh = invV*rho_ext; rhoh(2:Np,:)=0;
% rhoa = V*rhoh; rhoc = rhoa(1,:);
rhoc = Mesh.AVG1D*rho_ext;

% mmth = invV*mmt_ext; mmth(2:Np,:)=0;
% mmta = V*mmth; mmtc = mmta(1,:);
mmtc = Mesh.AVG1D*mmt_ext;

% Enerh = invV*Ener_ext; Enerh(2:Np,:)=0;
% Enera = V*Enerh; Enerc = Enera(1,:);
Enerc = Mesh.AVG1D*Ener_ext;

% Compute characterisic variables
if(~isempty(ind))
    for i=1:length(ind)
        
        [L,invL] = EulerCharMat1D(rhoc(1,ind(i)+1),mmtc(1,ind(i)+1),...
                                Enerc(1,ind(i)+1),gas_gamma,gas_const);
        
        Char  = invL*[rho_ext(:,ind(i))', rho_ext(:,ind(i)+1)', rho_ext(:,ind(i)+2)';
                      mmt_ext(:,ind(i))', mmt_ext(:,ind(i)+1)', mmt_ext(:,ind(i)+2)'
                      Ener_ext(:,ind(i))', Ener_ext(:,ind(i)+1)', Ener_ext(:,ind(i)+2)'];
        
        Char1 = SlopeLimit3(reshape(Char(1,:,:),[Mesh.Np,3]),Mesh.x(:,ind(i)),Limiter,Mesh);
        Char2 = SlopeLimit3(reshape(Char(2,:,:),[Mesh.Np,3]),Mesh.x(:,ind(i)),Limiter,Mesh);
        Char3 = SlopeLimit3(reshape(Char(3,:,:),[Mesh.Np,3]),Mesh.x(:,ind(i)),Limiter,Mesh);
        
        U     = L*[Char1(:)';Char2(:)'; Char3(:)' ];

        q_lim(:,ind(i),1) = reshape(U(1,:,:),[Mesh.Np,1]);
        q_lim(:,ind(i),2) = reshape(U(2,:,:),[Mesh.Np,1]);
        q_lim(:,ind(i),3) = reshape(U(3,:,:),[Mesh.Np,1]);        
    end
end


return;
