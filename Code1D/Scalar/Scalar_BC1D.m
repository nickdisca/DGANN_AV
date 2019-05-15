function [uin,uout] = Scalar_BC1D(bc_type,uL,uR)

if (strcmp(bc_type,'Periodic'))
    uin  = uR;
    uout = uL;
elseif(strcmp(bc_type,'Open'))
    uin  = uL;
    uout = uR;
else
    error('bc_type %s not available!!',bc_type)
end


return

