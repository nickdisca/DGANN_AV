function [Qf,neg_cells] = Euler_positivity_fix2D(Q,gas_gamma,AVG2D)

Qf = Q;
pre = Euler_Pressure2D(Q,gas_gamma);
neg_pcells = find(min(pre)<=0);
neg_dcells = find(min(Q(:,:,1))<=0);
neg_cells  = unique([neg_pcells,neg_dcells]);
for i=1:4
    Qf(:,neg_cells,i) = ones(length(Q(:,1)),1)*AVG2D*Q(:,neg_cells,i);
end

return