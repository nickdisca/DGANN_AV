function ulimit = SlopeLimit3(u,xl,Limiter,Mesh)

% Purpose: Apply slopelimiter (Pi^N) to u assuming u an N'th order polynomial 
% The input is the date for three cells only, with the limiting performed
% for the central cell. This funciton is used when the limiting needs
% to be done cell-by-cell

%Globals1D_DG;

[m,n] = size(u);
assert(n == 3);
ulimit = u;

if strcmp(Limiter,'NONE')
    return
else
    % Getting modal values
    uh = Mesh.invV*u; 
        
    % getting cell averages
    uavg = uh; uavg(2:end,:)=0; uavg = Mesh.V*uavg; v = uavg(1,:);
        
    % create piecewise linear solution for limiting on specified elements
    uh(3:end,:)=0; ul = Mesh.V*uh(:,2);
        
    % apply slope limiter to selected elements
    if strcmp(Limiter,'MINMOD')
        ulimit = SlopeLimitLin(ul,xl,v(1),v(2),v(3),Mesh);
    else
        error('Limiter %s not available!!',Limiter)
    end

end

return
