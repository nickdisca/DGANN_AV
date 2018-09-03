function Energy = Euler_Energy2D(rho,u,v,pre,gas_gamma)

Energy = 0.5*rho.*(u.^2 + v.^2) + pre/(gas_gamma-1);

return