function mu=Scalar2D_smooth_viscosity(mu_piece,Mesh)


mu_avg=1./Mesh.neighbors.*(Mesh.VToE*mu_piece'); 

mu=reshape(Mesh.Interp_matrix*mu_avg(Mesh.EToVT(:)),Mesh.Np,Mesh.K);

end