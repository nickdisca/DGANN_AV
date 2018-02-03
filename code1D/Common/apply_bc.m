function u_ext = apply_bc(u)

Globals1D_DG;

[m,n] = size(u);
u_ext = zeros(m,n+2);
u_ext(:,2:n+1) = u;
switch bc_type
    case 'Periodic'
        u_ext(:,1)   = u(:,n);
        u_ext(:,n+2) = u(:,1);
    case 'Open'
        u_ext(:,1)   = flipud(u(:,1));
        u_ext(:,n+2) = flipud(u(:,n));
    otherwise
        error('Uknown boundary condition %s',bc_type);
end

return