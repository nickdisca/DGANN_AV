function [entropy,ent_flux] = Set_scalar_entropy1D(model)

switch model
    case 'Advection'
        c = 1.0;
        entropy =@(u) 0.5*u.^2;
        ent_flux =@(u) c*0.5*u.^2;
    case 'Burgers'   
        entropy =@(u) 0.5*u.^2;
        ent_flux =@(u) 1/3*u.^3;
    case 'BuckLev'
        a = 0.5;
        entropy  =@(u) 0.5*u.^2;
        ent_flux =@(u) -a/(a+1)^2 * (  (a*(2-3*u)+u)./(a*(u-1).^2+u.^2) + log(u.^2 + a*(1-u).^2) +(a-1)/sqrt(a)*atan((a*(u-1)+u)/sqrt(a)) );  
    otherwise
        error('Unknown scalar model %s',model)
end

return