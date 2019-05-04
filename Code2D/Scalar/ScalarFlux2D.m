function [F,G] = ScalarFlux2D(Q,Problem)

switch Problem.model
    case 'Advection'
        F = Problem.AdvectionVelocity(1)*Q;
        G = Problem.AdvectionVelocity(2)*Q;
    case 'Burgers'
        F = Q.*Q/2;
        G = Q.*Q/2;    
    otherwise
        error('Unknown scalar model %s',Problem.model)
end

return
