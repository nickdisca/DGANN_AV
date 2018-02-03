function ind = Find_Tcells(u)

% Purpose: find all the troubled-cells for variable u
% NOTE: u must include ghost cell values

Globals1D_DG;
Globals1D_MLP;

% Extend based on bc_type
u_ext = apply_bc(u);

% Compute cell averages
uh = invV*u_ext; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);

eps0=1.0e-8;

% find end values of each element
ue1 = u(1,:); ue2 = u(end,:);

% find cell averages 
vk = v(2:K+1); vkm1 = v(1:K); vkp1 = v(3:K+2);

% Find elements in need of limiting
if(strcmp(indicator_type,'none'))
    ind = [];
elseif(strcmp(indicator_type,'minmod'))
    ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
    ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
    ind = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
elseif(strcmp(indicator_type,'TVB'))
    ve1 = vk - minmodB([(vk-ue1);vk-vkm1;vkp1-vk],TVB_M,x(end,:)-x(1,:));
    ve2 = vk + minmodB([(ue2-vk);vk-vkm1;vkp1-vk],TVB_M,x(end,:)-x(1,:));
    ind = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
elseif(strcmp(indicator_type,'NN'))
    bc_prob = ind_MLP([vkm1;vk;vkp1;ue1;ue2],n_input,n_output,...
        n_hidden_layer,leaky_alpha,WEIGHTS,BIASES);
    ind = find(bc_prob > 0.5-eps0);
else
    error('Indicator %s not available!!',indicator_type)
end


return
