function Q = IC(x, y)
 
% function Q = IC(x, y)
% Purpose: Set Initial conditio for 2D Advection. Simple sinecos wave

Q(:,:,1) = sin(4*pi*x).*cos(4*pi*y);
%Q(:,:,1) = 1+sin(2*pi*x).*sin(2*pi*y);
return;


