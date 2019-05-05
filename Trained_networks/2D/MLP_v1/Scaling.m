function x_in = Scaling(x)

% Scale each row my the maximum abs element of that row, provided the max
% abs is greater than 1. Else leave as is

[m,n]    = size(x);
fact     = repmat(max(max(abs(x)),1),m,1);
x_in     = x./fact;


return