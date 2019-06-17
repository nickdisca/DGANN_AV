function Visc_write1D(fid,t,mu)

fprintf(fid, '%.6f', t);
for k=1:length(mu(:))
    fprintf(fid, ', %.6f', mu(mod(k-1,size(mu,1))+1,floor((k-1)/size(mu,1))+1) );
end
fprintf(fid, '\n');

return