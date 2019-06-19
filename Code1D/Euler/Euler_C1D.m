function fval = Euler_C1D(u,v)

% Evauate local centered flux for 1D Euler equations
% u and v are conserved variables. It's just a simple arithmetic mean, but
% it's done for all components. It's NOT the centered flux for Euler
% equations!

fval=0.5*(u+v);

return