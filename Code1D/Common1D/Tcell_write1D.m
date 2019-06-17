function Tcell_write1D(fid,t,xc)

fprintf(fid, '%.6f', t);
for k=1:length(xc)
    fprintf(fid, ', %.6f', xc(k));
end
fprintf(fid, '\n');

return
