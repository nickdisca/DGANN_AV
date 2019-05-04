% Purpose : Compute inscribed circle diameter as characteristic
%           for grid to choose timestep

% Find vertex nodes
vmask1   = find( abs(Mesh.s+Mesh.r+2) < Mesh.NODETOL)'; 
vmask2   = find( abs(Mesh.r-1)   < Mesh.NODETOL)';
vmask3   = find( abs(Mesh.s-1)   < Mesh.NODETOL)';
vmask    = [vmask1;vmask2;vmask3]';
vx       = Mesh.x(vmask(:), :); vy = Mesh.y(vmask(:), :);

% Compute semi-perimeter and area
len1 = sqrt((vx(1,:)-vx(2,:)).^2+(vy(1,:)-vy(2,:)).^2);
len2 = sqrt((vx(2,:)-vx(3,:)).^2+(vy(2,:)-vy(3,:)).^2);
len3 = sqrt((vx(3,:)-vx(1,:)).^2+(vy(3,:)-vy(1,:)).^2);
sper = (len1 + len2 + len3)/2.0; 
Area = sqrt(sper.*(sper-len1).*(sper-len2).*(sper-len3));

% Compute scale using radius of inscribed circle
Mesh.dx = Area./sper;

% Compute scale using radius of circumscribed circle
Mesh.dx2 = len1.*len2.*len3./Area/4;