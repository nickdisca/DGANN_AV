function [re,se] = xy_ext_to_rs(x1,y1,x2,y2,x3,y3,xe,ye)

% function [re,se] = xy_ext_to_rs(Vcoor,xe,ye)
% Purpose : Changes xe,ye to re,se coordinates based on the triangle to
% right reference triangle transformation corresponding to the triangle 
% with vertex coordinates (x1,y1),(x2,y2),(x3,y3)

A = [-x1+x2 , -x1+x3;
     -y1+y2 , -y1+y3];

R = A\[2*xe' - x2-x3; 2*ye'-y2-y3];
re = R(1,:)';
se = R(2,:)';

return;
