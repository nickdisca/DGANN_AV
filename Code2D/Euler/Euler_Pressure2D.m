function Pressure = Euler_Pressure2D(Q,gas_gamma)

Pressure = (Q(:,:,4) - 0.5*(Q(:,:,2).^2 + Q(:,:,3).^2)./Q(:,:,1))*(gas_gamma-1);

return