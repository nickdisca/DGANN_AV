function Qlim = SlopeLimiter2D(Q,ind,Limiter)

switch Limiter
    
    case 'NONE'
        Qlim = Q;
    
    case 'BJES'
        Qlim = BarthJespLimiter2D(Q,ind);
        
    case 'VENK'
        Qlim = VenkLimiter2D(Q,ind);
        
    otherwise
        error('Unknown limiter option %s',Limiter)
end

return