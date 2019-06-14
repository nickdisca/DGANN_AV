function mu=Scalar1D_smooth_viscosity(mu_piece,x)

h=x(end,:)-x(1,:);
m=size(x,1)-1;

% piecewise linear smoothing
mu_piece_ext=[mu_piece(1), mu_piece, mu_piece(end)];
val_left=(mu_piece_ext(1:end-2)+mu_piece_ext(2:end-1))/2;
val_right=(mu_piece_ext(2:end-1)+mu_piece_ext(3:end))/2;
mu=ones(m+1,1)*val_left+(x-ones(m+1,1)*x(1,:))./h.*(ones(m+1,1)*(val_right-val_left));

end