function fval = ScalarLF1D(u,v,flux,dflux)
% Lax Friedrich dlux for scalar conservation laws

fval = 0.5*(flux(u) + flux(v) - max(abs(dflux(u)),abs(dflux(v))).*(v-u));

return
