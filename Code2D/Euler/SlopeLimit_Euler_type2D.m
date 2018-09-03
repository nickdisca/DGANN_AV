function [Q] = SlopeLimit_Euler_type2D(Q,gas_gamma,gas_const,Limiter,lim_var,ind)

% Limit solution based on type of limiting variable        

switch lim_var
    case 'prim'
        u        = Q(:,:,2)./Q(:,:,1);
        v        = Q(:,:,3)./Q(:,:,1);
        pre      = Euler_Pressure2D(Q,gas_gamma);
        
        Q(:,:,1) = SlopeLimiter2D(Q(:,:,1),ind,Limiter);
        u        = SlopeLimiter2D(u,ind,Limiter);
        v        = SlopeLimiter2D(v,ind,Limiter);
        pre      = SlopeLimiter2D(pre,ind,Limiter);
        
        Q(:,:,2) = Q(:,:,1).*u;
        Q(:,:,3) = Q(:,:,1).*v;
        Q(:,:,4) = Euler_Energy2D(Q(:,:,1),u,v,pre,gas_gamma);
    case 'con'
        Q(:,:,1) = SlopeLimiter2D(Q(:,:,1),ind,Limiter);
        Q(:,:,2) = SlopeLimiter2D(Q(:,:,2),ind,Limiter);
        Q(:,:,3) = SlopeLimiter2D(Q(:,:,3),ind,Limiter);
        Q(:,:,4) = SlopeLimiter2D(Q(:,:,4),ind,Limiter);
%     case 'char_cell'
%         [qc,Char_ext]    = EulerUtoC1D(q_ext,gas_gamma,gas_const);
%         Char(:,:,1)  = SlopeLimit1D(Char_ext(:,:,1),ind);
%         Char(:,:,2)  = SlopeLimit1D(Char_ext(:,:,2),ind);
%         Char(:,:,3)  = SlopeLimit1D(Char_ext(:,:,3),ind);
%         Q            = EulerCtoU1D(Char(:,:,:),qc(2:end-1,:),gas_gamma,gas_const);
%     case 'char_stencil'
%         Q = Euler_cstencil1D(q_ext,gas_gamma,gas_const,ind);
    otherwise
        error('Unknown limiting variable %s',lim_var) 
end

return
