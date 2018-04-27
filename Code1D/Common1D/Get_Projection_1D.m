% Purspose: Build Projection maps (needed for Shu-Fu indicator)
[m,n] = size(x);
x_ext = zeros(m,n+2);
x_ext(:,2:n+1) = x;

% Extend mesh depending on boundary condition
if(bc_cond{1,1} == 'P' || bc_cond{1,3} == 'P')
    x_ext(:,1)   = x(:,n)-(bnd_r-bnd_l);
    x_ext(:,n+2) = x(:,1)+(bnd_r-bnd_l);
else
    x_ext(:,1)   = x(:,1)-(x(end,1) - x(1,1));
    x_ext(:,n+2) = x(:,n)+(x(end,n) - x(1,n));
end

ProjectFromLeft1D  = zeros(Np,Np,K);
ProjectFromRight1D = zeros(Np,Np,K);
for i=1:K
    xl = (2*x_ext(:,i+1) - x_ext(1,i) -x_ext(end,i))/(x_ext(end,i)-x_ext(1,i));
    xr = (2*x_ext(:,i+1) - x_ext(1,i+2) -x_ext(end,i+2))/(x_ext(end,i+2)-x_ext(1,i+2));
    Vl = Vandermonde1D(N,xl); 
    Vr = Vandermonde1D(N,xr);
    ProjectFromLeft1D(:,:,i)  = Vl*invV;
    ProjectFromRight1D(:,:,i) = Vr*invV;
end
