function [q] = SlopeLimit_SWE_type(q,gravity,lim_var,ind)

% Limit solution based on type of limiting variable        

switch lim_var
    case 'prim'
        q(:,:,1) = SlopeLimit(q(:,:,1),ind);
        vel      = SlopeLimit(q(:,:,2)./q(:,:,1),ind);
        q(:,:,2) = q(:,:,1).*vel;
    case 'con'
        q(:,:,1) = SlopeLimit(q(:,:,1),ind);
        q(:,:,2) = SlopeLimit(q(:,:,2),ind);
    case 'char_cell'
        [qc,Char]    = SWEUtoC1D(q,gravity);
        Char(:,:,1)  = SlopeLimit(Char(:,:,1),ind);
        Char(:,:,2)  = SlopeLimit(Char(:,:,2),ind);
        q            = SWECtoU1D(Char,qc,gravity);
    case 'char_stencil'
        q = SWE_cstencil(q,gravity,ind);
    otherwise
        error('Unknown limiting variable %s',lim_var) 
end

return
