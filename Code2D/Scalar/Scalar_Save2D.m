% Saving plots to files

if(Save_soln)
   
    % Make solution dirctories if not existing
    mkdir('OUTPUT')
    
    fname = sprintf('OUTPUT/%s_DATA.mat',data_fname);
    save(fname,'x','y','invV','Save_times','Q_save','ind_save','ptc_hist','t_hist','sim_time')
    
   
end