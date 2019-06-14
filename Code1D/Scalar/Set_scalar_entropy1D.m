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



       if (type=='linear')
            E_old=1/2*uold.^2; E_new=1/2*u.^2;
            F_old=1/2*uold.^2; F_new=1/2*u.^2;
            f_prime_old=ones(size(uold)); f_prime_new=ones(size(u));
        end
        if (type=='burgers')
            E_old=1/2*uold.^2; E_new=1/2*u.^2;
            F_old=1/3*uold.^3; F_new=1/3*u.^3;
            f_prime_old=uold; f_prime_new=u;
        end
        if (type=='fourth')
            E_old=1/2*uold.^2; E_new=1/2*u.^2;
            F_old=1/5*uold.^5; F_new=1/5*u.^5;
            f_prime_old=uold.^3; f_prime_new=u.^3;
        end
        if (type=='buckley')
            E_old=1/2*uold.^2; E_new=1/2*u.^2;
            F_old=1/9*(2*(uold-2)./(3*uold.^2-2*uold+1)-2*log(3*uold.^2-2*uold+1)+sqrt(2)*atan((3*uold-1)/sqrt(2)));
            F_new=1/9*(2*(u-2)./(3*u.^2-2*u+1)-2*log(3*u.^2-2*u+1)+sqrt(2)*atan((3*u-1)/sqrt(2)));
            f_prime_old=4/9*uold.*(1-uold)./((uold.^2-2/3*uold+1/3).^2); 
            f_prime_new=4/9*u.*(1-u)./((u.^2-2/3*u+1/3).^2);
        end