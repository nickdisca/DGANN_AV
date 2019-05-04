function [TRI,xout,yout,interp] = GenInterpolators2D(Nout, Nin, xin, yin, invV)

% function [TRI,xout,yout,interp] = GenInterpolators2D(Nout, Nin, xin, yin)
% Purpose: Generate grid and interpolators for plotting fields
     
% build equally spaced grid on reference triangle (Right triangle)
Npout = (Nout+1)*(Nout+2)/2;
rout = zeros(Npout,1); sout = zeros(Npout,1); 
sk = 1;
for n=1:Nout+1
  for m=1:Nout+2-n
    rout(sk) = -1 + 2*(m-1)/Nout;
    sout(sk) = -1 + 2*(n-1)/Nout;
    counter(n,m) = sk; sk = sk+1;
  end
end

% build matrix to interpolate field data to equally spaced nodes
interp = InterpMatrix2D(rout, sout, Nin, invV);

% build triangulation of equally spaced nodes on reference triangle
tri = []; 
for n=1:Nout+1
  for m=1:Nout+1-n
    v1 = counter(n,m);   v2 = counter(n,m+1); 
    v3 = counter(n+1,m); v4 = counter(n+1,m+1);
    if(v4) 
      tri = [tri;[[v1 v2 v3];[v2 v4 v3]]]; 
    else
      tri = [tri;[[v1 v2 v3]]]; 
    end
  end
end

% build triangulation for all equally spaced nodes on all elements
% TRI = [];
% [dim1,dim2] = size(xin);
% for k=1:dim2
%   TRI = [TRI; tri+(k-1)*Npout];
% end
[dim1,dim2] = size(xin);
TRI   = repmat(tri,dim2,1);
shift = ones(length(tri(:,1)),1)*(Npout*(0:dim2-1));
TRI   = TRI + repmat(shift(:),1,3);


% interpolate node coordinates to equally spaced nodes
xout = interp*xin; yout = interp*yin;

return
