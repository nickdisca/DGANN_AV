function fval = ScalarC1D(u,v,flux,dflux)
% Centered flux for scalar conservation laws

fval = 0.5*(flux(u) + flux(v));

return