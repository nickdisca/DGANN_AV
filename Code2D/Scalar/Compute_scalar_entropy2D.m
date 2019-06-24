function [entropy,ent_flux_x,ent_flux_y] = Compute_scalar_entropy2D(Q,Problem)

switch Problem.model
    case 'Advection'
        c = Problem.AdvectionVelocity;
        entropy = Q.^2/2;
        ent_flux_x = c(1).*entropy;
        ent_flux_y = c(2).*entropy;
    case 'Burgers'
        entropy = Q.^2/2;
        ent_flux_x = Q.^3/3;
        ent_flux_y = Q.^3/3;        
    case 'Kpp'
        entropy = Q.^2/2;
        ent_flux_x = Q.*sin(Q)+cos(Q);
        ent_flux_y = Q.*cos(Q)-sin(Q);  
    otherwise
        error('Unknown scalar model %s',Problem.model)
end

return