% set up rendering
function Genplots_Euler2D(datapath_base,N,xran,yran,Plot_var,clines,viewval,tol)
close all
datapath = sprintf('%s_DATA',datapath_base);
load(datapath);

K  = length(x(1,:));
Np = length(x(:,1));
assert(Np==(N+1)*(N+2)/2,'Incorrect order prescribed for chosen data')

assert((iscell(clines) & length(clines) == length(Plot_var)),...
    'ERROR: ''clines'' must be a numeric cell of size equal to that of Plot_var')

%[TRI,xout,yout,interp] = GenInterpolators2D(N, N, x, y, invV);

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
plot(t_hist,ptc_hist,'LineWidth',2)
xlim([0,FinalTime])
ylim([0,100])
set(gca,'FontSize',15)
max_mark = max(ptc_hist);
avg_mark = mean(ptc_hist);
title(['Max cells marked = ',num2str(max_mark),'%, ',...
    'Avg cells marked = ',num2str(avg_mark),'%, ']);
% fname = sprintf('%s_tcell_perc_hist.pdf',datapath_base);
% print(fname,'-dpdf')
fname = sprintf('%s_tcell_perc_hist.png',datapath_base);
print(fname,'-dpng')

% figure(100)
% plot(t_hist,pnc_hist,'LineWidth',2)
% xlim([0,FinalTime])
% ylim([0,100])
% set(gca,'FontSize',20)
% max_mark = max(pnc_hist);
% avg_mark = mean(pnc_hist);
% [max_mark,avg_mark]*K/100
% K
% title(['Max cells marked = ',num2str(max_mark),'%, ',...
%     'Avg cells marked = ',num2str(avg_mark),'%, ']);
% % fname = sprintf('%s_tcell_perc_hist.pdf',datapath_base);
% % print(fname,'-dpdf')
% fname = sprintf('%s_ncell_perc_hist.png',datapath_base);
% print(fname,'-dpng')


nvar = length(Plot_var);



for t = 1:Nt
    
    figure(3)
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
    
    figure(4)
    plot(xavg(ind_save{1,t}), yavg(ind_save{1,t}), 'k.');
    axis equal
    xlim(xran)
    ylim(yran)
    xlabel('x')
    ylabel('y')
    %title(['Troubled cells at t=',num2str(FinalTime)]);
    set(gca,'FontSize',20)
    fname = sprintf('%s_tcell_t%d.pdf',datapath_base,t-1);
    print(fname,'-dpdf')
    
    for vind = 1:nvar
        
        vname = Plot_var{vind};
        
        QC = Q_save{1,t};
        
        switch vname
            case 'density'
                QC = QC(:,:,1);
            case 'velx'
                QC = QC(:,:,2)./QC(:,:,1);
            case 'vely'
                QC = QC(:,:,3)./QC(:,:,1);
            case 'pressure'
                QC = Euler_Pressure2D(QC,gas_gamma);
            case 'Energy'
                QC = QC(:,:,4);
            otherwise
                error(['Unknown variable type ',vname])
        end
        
        max(max(QC))
        min(min(QC))
        
        figure(4 + 2*vind-1)
        hold all
        F=scatteredInterpolant(x(:),y(:),QC(:)); interpc=F(xq,yq);
        [~,h] = contour(xq,yq,interpc,clines{1,vind});
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
        set(gca,'FontSize',15)
        
        figure(4 + 2*vind)
        %[xnew,ynew] = meshgrid(linspace(xran(1),xran(2),100),linspace(yran(1),yran(2),100));
        %interpc_new = interp2(xq,yq,interpc,xnew,ynew,'spline');
        %contour(xnew,ynew,interpc_new,clines{1,vind},'k-');
        contour(xq,yq,interpc,clines{1,vind},'k-');
        colormap('jet')
        axis equal
        xlim(xran)
        ylim(yran)
        xlabel('x')
        ylabel('y')
        set(gca,'FontSize',20)
        fname = sprintf('%s_%s_contour_bnw_t%d.pdf',datapath_base,vname,t-1);
        print(fname,'-dpdf')
%         fname = sprintf('%s_%s_contour_bnw_t%d.png',datapath_base,vname,t-1);
%         print(fname,'-dpng')
        
        
        
    end
end

figure(3)
grid on
fname = sprintf('%s_tcell_hist.pdf',datapath_base);
print(fname,'-dpdf')

for vind = 1:nvar
    
    vname = Plot_var{vind};
    
    QC = Q_save{1,t};
    
    switch vname
        case 'density'
            QC = QC(:,:,1);
        case 'velx'
            QC = QC(:,:,2)./QC(:,:,1);
        case 'vely'
            QC = QC(:,:,3)./QC(:,:,1);
        case 'pressure'
            QC = Euler_Pressure2D(QC,gas_gamma);
        case 'Energy'
            QC = QC(:,:,4);
        otherwise
            error(['Unknown variable type ',vname])
    end
    
    figure(4 + 2*vind-1)
    grid on
    fname = sprintf('%s_%s_contour_hist.pdf',datapath_base,vname);
    print(fname,'-dpdf')
    
%     figure(3 + nvar + vind)
%     PlotField2D(QC,interp,TRI,xout,yout); axis tight; drawnow;
%     xlim(xran)
%     ylim(yran)
%     xlabel('x')
%     ylabel('y')
%     zlabel('Solution')
%     set(gca,'FontSize',15)
    % ax = gca;
    % outerpos = ax.OuterPosition;
    % ti = ax.TightInset;
    % left = outerpos(1) + ti(1);
    % bottom = outerpos(2) + ti(2);
    % ax_width = outerpos(3) - ti(1) - ti(3);
    % ax_height = outerpos(4) - ti(2) - ti(4);
    % ax.Position = [left bottom ax_width ax_height];
    % fname = sprintf('OUTPUT/JPEG/%s_soln_3D.jpeg',data_fname);
    % print(fname,'-djpeg')
    
    
end



