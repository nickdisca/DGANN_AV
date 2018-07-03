function y = leaky_ReLU(x,alpha)

y = max(x,0.0) - alpha*max(-x,0.0);

return