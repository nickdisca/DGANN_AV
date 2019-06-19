function [entropy,ent_flux] = Set_Euler_entropy1D(gas_gamma)

pre = @(q) (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
                
entropy=@(q) q(:,:,1)/(gas_gamma-1).*log(pre(q)./(q(:,:,1).^gas_gamma));

ent_flux=@(q) q(:,:,2)./q(:,:,1).*entropy(q);

return