function Qlim = SlopeLimiter2D(Q,QG,ind,Limiter)

switch Limiter
    
    case 'NONE'
        Qlim = Q;
    
    case 'BJES'
        Qlim = BarthJespLimiter2D(Q,QG,ind);
        
    case 'VENK'
        Qlim = VenkLimiter2D(Q,QG,ind);
        
    otherwise
        error('Unknown limiter option %s',Limiter)
end

return