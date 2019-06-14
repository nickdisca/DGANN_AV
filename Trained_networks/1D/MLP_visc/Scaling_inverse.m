function out_scale = Scaling_inverse(out,u,h,wave_speed)

% Perform scaling for each column of out according to the physics of the
% problem, mulitplying by a mesh-dependent factor and a wave speed. Note
% that the output is a scalar value obtained via maximization.

[m,n] = size(out);
out=max(out);

u_ext = Apply_BC1D(u,bc_cond);
jump_L=u_ext(2,1:n)-u_ext(1,2:n+1);
jump_R=u_ext(2,2:n+1)-u_ext(1,3:n+2);
scaling=max(abs(jump_L),abs(jump_R));
scaling=min(scaling,h);
        
fact=scaling.*wave_speed; fact=repmat(fact,m,1);
out_scale = out.*fact;

return