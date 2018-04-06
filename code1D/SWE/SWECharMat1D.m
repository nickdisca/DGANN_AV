function [S,invS] = SWECharMat1D(depth,discharge,gravity)

S    = zeros(2,2);
invS = zeros(2,2);

S(1,1) = 1;
S(1,2) = 1;
S(2,1) = discharge./depth - sqrt(gravity*depth);
S(2,2) = discharge./depth + sqrt(gravity*depth);


Sdet      = 0.5./sqrt(gravity*depth);
% if(~isreal(Sdet))
%     fprintf('PASUINF');
%     pause(1000);
% end
invS(1,1) = (discharge./depth + sqrt(gravity*depth)).*Sdet;
invS(1,2) = -1.*Sdet;
invS(2,1) = (-discharge./depth + sqrt(gravity*depth)).*Sdet;
invS(2,2) = 1.*Sdet;


end
