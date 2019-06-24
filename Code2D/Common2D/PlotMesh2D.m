function PlotMesh2D(Mesh)

% function PlotMesh2D()
% Purpose: Show unstructured finite element grid

axis equal
xmax = max(max(Mesh.x)); xmin = min(min(Mesh.x));
ymax = max(max(Mesh.y)); ymin = min(min(Mesh.y));

Lx = xmax-xmin;
Ly = ymax-ymin;
xmax = xmax+.1*Lx; xmin = xmin-.1*Lx;
ymax = ymax+.1*Ly; ymin = ymin-.1*Ly;

axis([xmin xmax ymin ymax])
drawnow; pause(.05);

oFx = reshape(Mesh.Fx, Mesh.Nfp, Mesh.Nfaces*Mesh.K); oFy = reshape(Mesh.Fy, Mesh.Nfp, Mesh.Nfaces*Mesh.K);

figure(2)
plot(oFx, oFy, 'k-')
axis equal
axis([xmin xmax ymin ymax])

drawnow; pause(.05);
return;
