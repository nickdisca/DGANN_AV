function lvel = Get_scalar_vel2D(Q,model,c)

lvel = zeros(size(Q),2);

switch model
    case 'Advection'
        lvel(:,:,1) = c(1)*ones(size(Q));
        lvel(:,:,2) = c(2)*ones(size(Q));
    case 'Burgers'
        lvel(:,:,1) = Q;
        lvel(:,:,2) = Q;    
    otherwise
        error('Unknown scalar model %s',model)
end

return
