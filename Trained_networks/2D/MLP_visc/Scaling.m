function x_in = Scaling(x)

% Perform scaling for each column of x in such a way that each sample is in
% [-1,1]. 

[m,n] = size(x);
fact=max(abs(x)); fact=repmat(fact,m,1)+1e-8;
x_in = x./fact;

return