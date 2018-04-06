function [flux,dflux] = Set_scalar_flux1D(model)

switch model
    case 'Advection'
        c = 1.0;
        flux =@(u) c*u;
        dflux =@(u) c*ones(size(u));
    case 'Burgers'   
        flux =@(u) 0.5*u.^2;
        dflux =@(u) u;
    case 'BuckLev'
        a = 0.5;
        flux  =@(u) u.^2./(u.^2 + a*(1-u).^2);
        dflux =@(u) 2*a*u.*(1-u)./(u.^2 + a*(1-u).^2).^2;  
    otherwise
        error('Unknown scalar model %s',model)
end

return
