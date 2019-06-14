function [y]=softplus(x)

y=log(1+exp(x));

indexes=(x>=30);
y(indexes)=x(indexes);

end