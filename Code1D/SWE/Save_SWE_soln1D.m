function Save_SWE_soln1D(q,x,fname_base)

% Globals1D_DG;
% Globals1D_MLP;

fid = fopen(strcat(fname_base,'.dat'),'w');                
depth = q(:,:,1);
vel   = q(:,:,2)./q(:,:,1);
fprintf(fid, '%.16f %.16f %.16f \n', [x(:) depth(:) vel(:)]');
fclose(fid);


return
