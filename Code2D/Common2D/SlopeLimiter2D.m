function Qlim = SlopeLimiter2D(Q,QG,ind,Limiter,Mesh)

switch Limiter
    
    case 'NONE'
        Qlim = Q;
    
    case 'BJES'
        Qlim = BarthJespLimiter2D(Q,QG,ind,Mesh);
        
    case 'VENK'
        Qlim = VenkLimiter2D(Q,QG,ind,Mesh);   
        
    otherwise
        error('Unknown limiter option %s',Limiter)
end

return