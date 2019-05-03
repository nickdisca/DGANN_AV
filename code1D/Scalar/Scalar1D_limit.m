function u = Scalar1D_limit(u,ind,bc_cond,Limiter,Mesh)

u_ext = Apply_BC1D(u,bc_cond);
u     = SlopeLimit1D(u_ext,ind,Limiter,Mesh);

return
