% Saving data to files

if(Output.save_soln)
    
    % Make solution dirctories if not existing
    mkdir('OUTPUT')
    
    fname = sprintf('OUTPUT/%s_DATA.mat',data_fname);
    x = Mesh.x;
    y = Mesh.y;
    invV = Mesh.invV;
    
    
    save(fname,'x','y','invV','Save_times','Q_save','ind_save','ptc_hist','pnc_hist','t_hist','sim_time','gas_gamma','gas_const')
    
    
end
