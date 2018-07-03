function x_in = Scaling(x)

% Perform scaling for each column of x. If the values in a column do not 
% lie in [-1,1], then the column is divided by the maximum of the absolute 
% values in that column.

[m,n] = size(x);
fact = repmat(max(max(abs(x),[],1),1),m,1);
x_in = x./fact;

return