function maxeig = Get_scalar_eig2D(Q,model,c)

switch model
    case 'Advection'
        maxeig = sqrt(c(1)^2 + c(2)^2)*ones(size(Q));
    case 'Burgers'
        maxeig = abs(Q);    
    otherwise
        error('Unknown scalar model %s',model)
end

return
