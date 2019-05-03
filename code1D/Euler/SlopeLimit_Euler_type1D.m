function [q] = SlopeLimit_Euler_type1D(q,ind,Problem,Limit,Mesh)

% Limit solution based on type of limiting variable        

% bc = {'N',0.0,'N',0.0};

q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),Problem.bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),Problem.bc_cond(2,:));
q_ext(:,:,3)  = Apply_BC1D(q(:,:,3),Problem.bc_cond(3,:));

switch Limit.lim_var
    case 'prim'
        vel_ext  = q_ext(:,:,2)./q_ext(:,:,1);
        pre_ext  = (Problem.gas_gamma-1)*(q_ext(:,:,3) - 0.5*q_ext(:,:,2).^2./q_ext(:,:,1));
        
        q(:,:,1) = SlopeLimit1D(q_ext(:,:,1),ind,Limit.Limiter,Mesh);
        vel      = SlopeLimit1D(vel_ext,ind,Limit.Limiter,Mesh);
        pre      = SlopeLimit1D(pre_ext,ind,Limit.Limiter,Mesh);
        
        q(:,:,2) = q(:,:,1).*vel;
        q(:,:,3) = 0.5*q(:,:,1).*vel.^2 + pre/(Problem.gas_gamma-1);
    case 'con'
        q(:,:,1) = SlopeLimit1D(q_ext(:,:,1),ind,Limit.Limiter,Mesh);
        q(:,:,2) = SlopeLimit1D(q_ext(:,:,2),ind,Limit.Limiter,Mesh);
        q(:,:,3) = SlopeLimit1D(q_ext(:,:,3),ind,Limit.Limiter,Mesh);
    case 'char_cell'
        [qc,Char_ext]    = EulerUtoC1D(q_ext,Problem.gas_gamma,Problem.gas_const,Mesh);
        Char(:,:,1)      = SlopeLimit1D(Char_ext(:,:,1),ind,Limit.Limiter,Mesh);
        Char(:,:,2)      = SlopeLimit1D(Char_ext(:,:,2),ind,Limit.Limiter,Mesh);
        Char(:,:,3)      = SlopeLimit1D(Char_ext(:,:,3),ind,Limit.Limiter,Mesh);
        q                = EulerCtoU1D(Char(:,:,:),qc(2:end-1,:),Problem.gas_gamma,Problem.gas_const,Mesh);
    case 'char_stencil'
        q = Euler_cstencil1D(q_ext,Problem.gas_gamma,Problem.gas_const,ind,Limit.Limiter,Mesh);
    otherwise
        error('Unknown limiting variable %s',Limit.lim_var) 
end

return
