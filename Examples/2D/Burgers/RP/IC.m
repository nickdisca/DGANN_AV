function Q = IC(x, y)
 
% function Q = IC(x, y)
% Purpose: Set Initial conditio for 2D Burgers.

Q(:,:,1) = -1*(x>0).*(y>0) + ...
           -0.2*(x<=0).*(y>0) + ...
           0.5*(x<=0).*(y<=0) + ...
           0.8*(x>0).*(y<=0);

% Q(:,:,1) = 1 * 100*x + 20*y;

return;


