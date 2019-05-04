function Q = IC(x, y)
 
% function Q = IC(x, y)
% Purpose: Set Initial conditio for 2D Burgers.

Q(:,:,1) = sin(2*pi*(x+0.5)).*sin(2*pi*(y+0.5)).*(abs(x)<=0.5).*(abs(y)<=0.5);

% Q(:,:,1) = -1*(x>0).*(y>0) + ...
%            -0.2*(x<=0).*(y>0) + ...
%            0.5*(x<=0).*(y<=0) + ...
%            0.8*(x>0).*(y<=0);
return;


