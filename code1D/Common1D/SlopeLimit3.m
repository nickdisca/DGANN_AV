function ulimit = SlopeLimit3(u,xl)

% Purpose: Apply slopelimiter (Pi^N) to u assuming u an N'th order polynomial 
% The input is the date for three cells only, with the limiting performed
% for the central cell. This funciton is used when the limiting needs
% to be done cell-by-cell

Globals1D_DG;

[m,n] = size(u);
assert(n == 3);
ulimit = u;

if strcmp(rec_limiter,'none')
    return
else
    % Getting modal values
    uh = invV*u; 
        
    % getting cell averages
    uavg = uh; uavg(2:Np,:)=0; uavg = V*uavg; v = uavg(1,:);
        
    % create piecewise linear solution for limiting on specified elements
    uh(3:Np,:)=0; ul = V*uh(:,2);
        
    % apply slope limiter to selected elements
    if strcmp(rec_limiter,'minmod')
        ulimit = SlopeLimitLin(ul,xl,v(1),v(2),v(3));
    else
        error('Limiter %s not available!!',rec_limiter)
    end

end

return
