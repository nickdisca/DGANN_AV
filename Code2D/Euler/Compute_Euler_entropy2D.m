function [entropy,ent_flux_X,ent_flux_Y] = Compute_Euler_entropy2D(Q,gas_gamma)

    vel_u = Q(:,:,2)./Q(:,:,1); vel_v = Q(:,:,3)./Q(:,:,1);
    pres = (gas_gamma-1.0)*(Q(:,:,4) - Q(:,:,1).*(vel_u.^2+vel_v.^2)/2); 
                
    entropy= Q(:,:,1)/(gas_gamma-1).*log(pres./(Q(:,:,1).^gas_gamma));

    ent_flux_X= Q(:,:,2)./Q(:,:,1).*entropy;
    ent_flux_Y= Q(:,:,3)./Q(:,:,1).*entropy;

return