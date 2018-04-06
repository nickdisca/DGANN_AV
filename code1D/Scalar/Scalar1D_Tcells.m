function ind = Scalar1D_Tcells(u,bc_cond)

u_ext = Apply_BC1D(u,bc_cond);
ind   = Find_Tcells1D(u_ext);

return

