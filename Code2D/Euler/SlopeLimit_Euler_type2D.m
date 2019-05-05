function [Q] = SlopeLimit_Euler_type2D(Q,QG,ind,gas_gamma,Limit,Mesh)

% Limit solution based on type of limiting variable        

switch Limit.lim_var
    case 'prim'
        u        = Q(:,:,2)./Q(:,:,1);
        v        = Q(:,:,3)./Q(:,:,1);
        pre      = Euler_Pressure2D(Q,gas_gamma);
        
        uG       = QG(:,:,2)./QG(:,:,1);
        vG       = QG(:,:,3)./QG(:,:,1);
        preG     = Euler_Pressure2D(QG,gas_gamma);
        
        Q(:,:,1) = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind,Limit.Limiter,Mesh);
        u        = SlopeLimiter2D(u,uG,ind,Limit.Limiter,Mesh);
        v        = SlopeLimiter2D(v,vG,ind,Limit.Limiter,Mesh);
        pre      = SlopeLimiter2D(pre,preG,ind,Limit.Limiter,Mesh);
        
        Q(:,:,2) = Q(:,:,1).*u;
        Q(:,:,3) = Q(:,:,1).*v;
        Q(:,:,4) = Euler_Energy2D(Q(:,:,1),u,v,pre,gas_gamma);
    case 'con'
        Q(:,:,1) = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind,Limit.Limiter,Mesh);
        Q(:,:,2) = SlopeLimiter2D(Q(:,:,2),QG(:,:,2),ind,Limit.Limiter,Mesh);
        Q(:,:,3) = SlopeLimiter2D(Q(:,:,3),QG(:,:,3),ind,Limit.Limiter,Mesh);
        Q(:,:,4) = SlopeLimiter2D(Q(:,:,4),QG(:,:,4),ind,Limit.Limiter,Mesh);
    otherwise
        error('Unknown limiting variable %s',Limit.lim_var) 
end


return
