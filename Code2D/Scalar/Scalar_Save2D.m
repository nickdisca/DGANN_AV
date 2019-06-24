% Saving plots to files

% Make solution dirctories if not existing
mkdir('OUTPUT')

fname = sprintf('OUTPUT/%s_DATA.mat',data_fname);
x = Mesh.x;
y = Mesh.y;
invV = Mesh.invV;
J = Mesh.J;
Mass = Mesh.MassMatrix;

save(fname,'x','y','invV','J','Mass','Save_times','Q_save','ind_save','visc_save','ptc_hist','maxvisc_hist','t_hist','sim_time')


