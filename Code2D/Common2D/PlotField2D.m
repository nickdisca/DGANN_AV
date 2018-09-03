function PlotField2D(uin, interp,TRI,xout,yout)

% function uout = PlotField2D(interp,TRI,xout,yout)
% Purpose: filled contour plot of solution data

% interpolate field to equally spaced nodes
uout = interp*uin;


% render and format solution field
trisurf(TRI, xout(:), yout(:), uout(:));
shading interp,    material shiny,    lighting gouraud 
camlight headlight
return
