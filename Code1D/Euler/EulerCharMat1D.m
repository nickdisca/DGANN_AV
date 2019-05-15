function [S,invS] = EulerCharMat1D(rho,rhou,Ener,gas_gamma,gas_const)

u = rhou./rho;
p = (Ener - 0.5*u.^2.*rho)*(gas_gamma-1);
a = sqrt(gas_gamma*p./rho);
H = (Ener + p)./rho;

S    = zeros(3,3);
invS = zeros(3,3);

S(1,1) = 1./(2*a);
S(1,2) = 1;
S(1,3) = 1./(2*a);
S(2,1) = (u + a)./(2*a);
S(2,2) = u;
S(2,3) = (u - a)./(2*a);
S(3,1) = (H + u.*a)./(2*a);
S(3,2) = u.^2/2;
S(3,3) = (H - u.*a)./(2*a);


invS(1,1) = (0.5*(gas_gamma - 1)*u.^2 - a.*u)./a;
invS(1,2) = -((gas_gamma - 1)*u - a)./a;
invS(1,3) = (gas_gamma - 1)./a;
invS(2,1) = 1 - 0.5*(gas_gamma - 1)*u.^2./a.^2;
invS(2,2) = (gas_gamma - 1)./a.^2.*u;
invS(2,3) = -(gas_gamma - 1)./a.^2;
invS(3,1) = (0.5*(gas_gamma - 1)*u.^2 + a.*u)./a;
invS(3,2) = -((gas_gamma - 1)*u + a)./a;
invS(3,3) = (gas_gamma - 1)./a;


end
