function Save_scalar_soln1D(u,x,fname_base)

% Globals1D_DG;
% Globals1D_MLP;

fprintf('... saving solution data files in directory OUTPUT\n')
fid = fopen(strcat(fname_base,'.dat'),'w');
fprintf(fid, '%.16f %.16f \n', [x(:) u(:)]');
fclose(fid);


return
