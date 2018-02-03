function [q] = SlopeLimit_Euler_type(q,gas_gamma,gas_const,lim_var,ind)

% Limit solution based on type of limiting variable        

switch lim_var
    case 'prim'
        vel      = q(:,:,2)./q(:,:,1);
        pre      = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
        
        q(:,:,1) = SlopeLimit(q(:,:,1),ind);
        vel      = SlopeLimit(vel,ind);
        pre      = SlopeLimit(pre,ind);
        
        q(:,:,2) = q(:,:,1).*vel;
        q(:,:,3) = 0.5*q(:,:,1).*vel.^2 + pre/(gas_gamma-1);
    case 'con'
        q(:,:,1) = SlopeLimit(q(:,:,1),ind);
        q(:,:,2) = SlopeLimit(q(:,:,2),ind);
        q(:,:,2) = SlopeLimit(q(:,:,2),ind);
    case 'char_cell'
        [qc,Char]    = EulerUtoC1D(q,gas_gamma,gas_const);
        Char(:,:,1)  = SlopeLimit(Char(:,:,1),ind);
        Char(:,:,2)  = SlopeLimit(Char(:,:,2),ind);
        Char(:,:,3)  = SlopeLimit(Char(:,:,3),ind);
        q            = EulerCtoU1D(Char,qc,gas_gamma,gas_const);
    case 'char_stencil'
        q = Euler_cstencil(q,gas_gamma,gas_const,ind);
    otherwise
        error('Unknown limiting variable %s',lim_var) 
end

return
