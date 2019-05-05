function PlotEulerFields2D(Q,time,gas_gamma,gas_const,x,y,xq,yq,TRI,xout,yout,interp,Output,figc,figf)

nvar  = numel(Output.plot_var);

if(nvar ==1)
    nrows = 1;
    ncols = 1;
else
    nrows = 2;
    ncols = ceil(nvar/2);
end

for n = 1:nvar
    if(strcmp(Output.plot_var{n},'density'))
        QC = Q(:,:,1);
    elseif(strcmp(Output.plot_var{n},'velx'))
        QC = Q(:,:,2)./Q(:,:,1);
    elseif(strcmp(Output.plot_var{n},'vely'))
        QC = Q(:,:,3)./Q(:,:,1);
    elseif(strcmp(Output.plot_var{n},'pressure'))
        QC = Euler_Pressure2D(Q,gas_gamma);  
    elseif(strcmp(Output.plot_var{n},'energy'))
        QC = Q(:,:,4);   
    end
    
    figure(figc)
    subplot(nrows,ncols,n)
    F  = scatteredInterpolant(x(:),y(:),QC(:)); interpc=F(xq,yq);
    contour(xq,yq,interpc,Output.clines{n},'k-','LineWidth',1);
    xlim(Output.xran)
    ylim(Output.yran)
    xlabel('x')
    ylabel('y')
    title([Output.plot_var{n},' at t=',num2str(time)])
    drawnow;
    

    figure(figf)
    subplot(nrows,ncols,n)
    PlotField2D(QC,interp,TRI,xout,yout);
    xlim(Output.xran)
    ylim(Output.yran)
    zlabel(Output.plot_var{n})
    xlabel('x')
    ylabel('y')
    title(['t=',num2str(time)])
    drawnow;
    
end