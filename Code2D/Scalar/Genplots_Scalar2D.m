% set up rendering
function Genplots_Scalar2D(datapath_base,N,xran,yran,clines,viewval,tol,TRI,xout,yout,interp)
close all
datapath = sprintf('%s_DATA',datapath_base);
load(datapath);

K  = length(x(1,:));
Np = length(x(:,1));
assert(Np==(N+1)*(N+2)/2,'Incorrect order prescribed for chosen data')

MassMatrix = invV'*invV;
AVG2D      = sum(MassMatrix)/2;
xavg = AVG2D*x;
yavg = AVG2D*y;

if(abs(Save_times(end-1)-Save_times(end))<tol)
    Save_times = [Save_times(1:end-2),Save_times(end)];
end

Nt        = length(Save_times);
FinalTime = Save_times(end);

% Plotting variables
Nq = 200;
hx = (xran(2)-xran(1))/Nq; hy = (yran(2)-yran(1))/Nq;
[xq,yq]=meshgrid(xran(1):hx:xran(2),yran(1):hy:yran(2));

zoom_ratio = 0.5*Nt/FinalTime*min(xran(2)-xran(1),yran(2)-yran(1));

figure(1)
xlim([0,FinalTime])
ylim([0,100])
plot(t_hist,ptc_hist,'LineWidth',2)
set(gca,'FontSize',20)
max_mark = max(ptc_hist);
avg_mark = mean(ptc_hist);
title(['Max cells marked = ',num2str(max_mark),'%, ',...
       'Avg cells marked = ',num2str(avg_mark),'%, ']);
fname = sprintf('%s_tcell_perc_hist.pdf',datapath_base);
print(fname,'-dpdf')


for t = 1:Nt
    
    figure(2)
    hold all
    xlim(xran)
    ylim(yran)
    xlabel('x')
    ylabel('y')
    zlabel('time')
    plot3(xavg(ind_save{1,t}), yavg(ind_save{1,t}),Save_times(t)*ones(1,length(ind_save{1,t})), 'k.');
    zticks(Save_times);
    set(gca,'FontSize',20)
    ax = gca;
    ax.DataAspectRatio = [zoom_ratio zoom_ratio 1];
    view(viewval)
    
    figure(3+2*t-1)
    plot(xavg(ind_save{1,t}), yavg(ind_save{1,t}), 'k.');
    xlim(xran)
    ylim(yran)
    xlabel('x')
    ylabel('y')
    %title(['Troubled cells at t=',num2str(FinalTime)]);
    set(gca,'FontSize',20)
    fname = sprintf('%s_tcell_t%d.pdf',datapath_base,t-1);
    print(fname,'-dpdf')
    
    % fname = sprintf('OUTPUT/%s_tcells.pdf',data_fname);
    % print(fname,'-dpdf')
    
    QC = Q_save{1,t};
    
    figure(3)
    hold all
    F=scatteredInterpolant(x(:),y(:),QC(:)); interpc=F(xq,yq);
    [~,h] = contour(xq,yq,interpc,clines);
    h.ContourZLevel = Save_times(t);
    xlim(xran)
    ylim(yran)
    zticks(Save_times);
    xlabel('x')
    ylabel('y')
    zlabel('time')
    colormap(gray)
    ax = gca;
    ax.DataAspectRatio = [zoom_ratio zoom_ratio 1];
    view(viewval)
    set(gca,'FontSize',20)
    
    figure(3+2*t)
    F=scatteredInterpolant(x(:),y(:),QC(:)); interpc=F(xq,yq);
    contour(xq,yq,interpc,clines,'k-');
    xlim(xran)
    ylim(yran)
    xlabel('x')
    ylabel('y')
    set(gca,'FontSize',20)
    fname = sprintf('%s_contour_bnw_t%d.pdf',datapath_base,t-1);
    print(fname,'-dpdf')
    
    

end

figure(2)
grid on
fname = sprintf('%s_tcell_hist.pdf',datapath_base);
print(fname,'-dpdf')

figure(3)
grid on
fname = sprintf('%s_contour_hist.pdf',datapath_base);
print(fname,'-dpdf')

figure(100)
grid on
PlotField2D(QC,interp,TRI,xout,yout); axis tight; drawnow;
fname = sprintf('%s_soln_final.pdf',datapath_base);
print(fname,'-dpdf')




