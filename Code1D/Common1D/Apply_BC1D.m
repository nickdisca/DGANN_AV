function u_ext = Apply_BC1D(u,bc)

[m,n] = size(u);
u_ext = zeros(m,n+2);
u_ext(:,2:n+1) = u;

% Periodic extension
if(bc{1,1} == 'P' || bc{1,3} == 'P')
    u_ext(:,1)   = u(:,n);
    u_ext(:,n+2) = u(:,1);
    return
end

% Left extension
if(bc{1,1} == 'N')
    u_ext(:,1)   = flipud(u(:,1));
elseif(bc{1,1} == 'D')
    u_ext(:,1)   = -flipud(u(:,1)) + 2*bc{1,2};
end

% Right extension
if(bc{1,3} == 'N')
    u_ext(:,n+2) = flipud(u(:,n));
elseif(bc{1,3} == 'D')
    u_ext(:,n+2)   = -flipud(u(:,n)) + 2*bc{1,4};
end


return
