function [IM] = InterpMatrix2D(rout, sout, N, invV)

% function [IM] = InterpMatrix2D(rout, sout, N, invV)
% Purpose: Compute local elemental interpolation matrix
 
% compute Vandermonde at (rout,sout)
Vout = Vandermonde2D(N, rout, sout);

% build interpolation matrix
IM = Vout*invV;
return

  
  
