function out_scale = Scaling_inverse_2D(out,u,h,wave_speed,Mesh)

% Perform scaling for each column of out according to the physics of the
% problem, mulitplying by a mesh-dependent factor and a wave speed. Note
% that the output is a scalar value obtained via maximization.

[m,n] = size(out);
out=max(out);

scaling=max(reshape(abs(u(Mesh.vmapP)-u(Mesh.vmapM)),3*(Mesh.N+1),Mesh.K));
scaling=min(scaling,h);

fact=scaling.*wave_speed;
out_scale = out.*fact;

return