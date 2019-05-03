function [q] = SlopeLimit_SWE_type1D(q,ind,gravity,bc_cond,Limit,Mesh)

% Limit solution based on type of limiting variable        

q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),bc_cond(2,:));

switch Limit.lim_var
    case 'prim'
        q(:,:,1) = SlopeLimit1D(q_ext(:,:,1),ind,Limit.Limiter,Mesh);
        vel      = SlopeLimit1D(q_ext(:,:,2)./q_ext(:,:,1),ind,Limit.Limiter,Mesh);
        q(:,:,2) = q(:,:,1).*vel;
    case 'con'
        q(:,:,1) = SlopeLimit1D(q_ext(:,:,1),ind,Limit.Limiter,Mesh);
        q(:,:,2) = SlopeLimit1D(q_ext(:,:,2),ind,Limit.Limiter,Mesh);
    case 'char_cell'
        [qc,Char]    = SWEUtoC1D(q_ext,gravity,Mesh);
        Char(:,:,1)  = SlopeLimit1D(Char(:,:,1),ind,Limit.Limiter,Mesh);
        Char(:,:,2)  = SlopeLimit1D(Char(:,:,2),ind,Limit.Limiter,Mesh);
        q            = SWECtoU1D(Char(:,2:Mesh.K-1,:),qc(2:Mesh.K+1,:),gravity);
    case 'char_stencil'
        q = SWE_cstencil1D(q_ext,gravity,ind,Limit.Limiter,Mesh);
    otherwise
        error('Unknown limiting variable %s',Limit.lim_var) 
end

return
