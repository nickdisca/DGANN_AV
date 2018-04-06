function [q] = SlopeLimit_Euler_type1D(q,gas_gamma,gas_const,lim_var,ind,bc_cond)

% Limit solution based on type of limiting variable        

% bc = {'N',0.0,'N',0.0};

q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),bc_cond(2,:));
q_ext(:,:,3)  = Apply_BC1D(q(:,:,3),bc_cond(3,:));

switch lim_var
    case 'prim'
        vel_ext  = q_ext(:,:,2)./q_ext(:,:,1);
        pre_ext  = (gas_gamma-1)*(q_ext(:,:,3) - 0.5*q_ext(:,:,2).^2./q_ext(:,:,1));
        
        q(:,:,1) = SlopeLimit1D(q_ext(:,:,1),ind);
        vel      = SlopeLimit1D(vel_ext,ind);
        pre      = SlopeLimit1D(pre_ext,ind);
        
        q(:,:,2) = q(:,:,1).*vel;
        q(:,:,3) = 0.5*q(:,:,1).*vel.^2 + pre/(gas_gamma-1);
    case 'con'
        q(:,:,1) = SlopeLimit1D(q_ext(:,:,1),ind);
        q(:,:,2) = SlopeLimit1D(q_ext(:,:,2),ind);
        q(:,:,3) = SlopeLimit1D(q_ext(:,:,3),ind);
    case 'char_cell'
        [qc,Char_ext]    = EulerUtoC1D(q_ext,gas_gamma,gas_const);
        Char(:,:,1)  = SlopeLimit1D(Char_ext(:,:,1),ind);
        Char(:,:,2)  = SlopeLimit1D(Char_ext(:,:,2),ind);
        Char(:,:,3)  = SlopeLimit1D(Char_ext(:,:,3),ind);
        q            = EulerCtoU1D(Char(:,:,:),qc(2:end-1,:),gas_gamma,gas_const);
    case 'char_stencil'
        q = Euler_cstencil1D(q_ext,gas_gamma,gas_const,ind);
    otherwise
        error('Unknown limiting variable %s',lim_var) 
end

return
