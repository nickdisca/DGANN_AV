function maxeig = Get_scalar_eig2D(Q,Problem)

switch Problem.model
    case 'Advection'
        c = Problem.AdvectionVelocity;
        maxeig = sqrt(c(1)^2 + c(2)^2)*ones(size(Q));
    case 'Burgers'
        maxeig = abs(Q);    
    otherwise
        error('Unknown scalar model %s',Problem.model)
end

return
