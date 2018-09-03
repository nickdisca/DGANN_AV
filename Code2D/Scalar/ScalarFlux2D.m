function [F,G] = ScalarFlux2D(Q,model,c)

switch model
    case 'Advection'
        F = c(1)*Q;
        G = c(2)*Q;
    case 'Burgers'
        F = Q.*Q/2;
        G = Q.*Q/2;    
    otherwise
        error('Unknown scalar model %s',model)
end

return
