function ind = Find_Tcells1D(u,Mesh,Limit,Net)

% Purpose: find all the troubled-cells for variable u
% NOTE: u must include ghost cell values

% Globals1D_DG;
% Globals1D_MLP;

% Compute cell averages
%uh = invV*u; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);
v = Mesh.AVG1D*u;

eps0=1.0e-8;

% find end values of each element (excluding ghost cells)
ue1 = u(1,2:end-1); ue2 = u(end,2:end-1);

% find cell averages 
vk = v(2:Mesh.K+1); vkm1 = v(1:Mesh.K); vkp1 = v(3:Mesh.K+2);

% Find elements in need of limiting
if(strcmp(Limit.Indicator,'NONE'))
    ind = [];
elseif(strcmp(Limit.Indicator,'ALL'))
    ind = 1:Mesh.K;    
elseif(strcmp(Limit.Indicator,'MINMOD'))
    ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
    ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
    ind = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
elseif(strcmp(Limit.Indicator,'TVB'))
    ve1 = vk - minmodB([(vk-ue1);vk-vkm1;vkp1-vk],Limit.TVBM,Mesh.x(end,:)-Mesh.x(1,:));
    ve2 = vk + minmodB([(ue2-vk);vk-vkm1;vkp1-vk],Limit.TVBM,Mesh.x(end,:)-Mesh.x(1,:));
    ind = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
elseif(strcmp(Limit.Indicator,'NN'))
    assert(Net.avail)
    bc_prob = ind_MLP1D([vkm1;vk;vkp1;ue1;ue2],Net);
    ind = find(bc_prob > 0.5-eps0);
% elseif(strcmp(Indicator,'FuShu'))
%     ind_val = zeros(1,K);
%     for i = 1:K % Need to replace loop by faster algorithm
%         upl = ProjectFromLeft1D(:,:,i)*u(:,i);
%         %uhl = invV*upl; uhl(2:end) = 0; ulavg = V*uhl; vl = ulavg(1);
%         vl  = AVG1D*upl;
%         upr = ProjectFromRight1D(:,:,i)*u(:,i+2);
%         %uhr = invV*upr; uhr(2:end) = 0; uravg = V*uhr; vr = uravg(1);
%         vr  = AVG1D*upr;
% %         figure(100)
% %         plot(x(:,i),u(:,i+1),'r--o')
% %         hold all
% %         plot([x(:,i)-(x(end,i)-x(1,i));x(:,i)],[u(:,i);upl],'b--x')
% %         plot([x(:,i);x(:,i)+(x(end,i)-x(1,i))],[upr;u(:,i+2)],'g--s')
% %         hold off
%         ind_val(i) = (abs(vk(i) - vl) + abs(vk(i) - vr))/max(abs([vk(i),vkm1(i),vkp1(i),eps]));
%     end
%     ind = find(ind_val > ck1D(N));  
else
    error('Indicator %s not available!!',Limit.Indicator)
end


return
