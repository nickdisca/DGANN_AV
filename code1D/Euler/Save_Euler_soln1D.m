function Save_Euler_soln1D(q,x,gas_gamma,fname_base)



fid = fopen(strcat(fname_base,'.dat'),'w');  
rho = q(:,:,1);
vel = q(:,:,2)./q(:,:,1);
pre = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
fprintf(fid, '%.16f %.16f %.16f  %.16f\n', [x(:) rho(:) vel(:) pre(:)]');
fclose(fid);


return
