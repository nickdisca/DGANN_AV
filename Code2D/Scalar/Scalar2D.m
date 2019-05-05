function [Q_save,ind_save,ptc_hist,t_hist,Save_times] = Scalar2D(Q,Problem,Mesh,Limit,Net,Output)

% function Q = Scalar2D(Q)
% Purpose  : Integrate 2D Scalar equation using a SSP-RK3


% Get centroid
xavg = Mesh.AVG2D*Mesh.x;
yavg = Mesh.AVG2D*Mesh.y;

% Get max wave speed
maxeig = Get_scalar_eig2D(Q,Problem);

% compute initial timestep
if(isempty(Problem.fixed_dt))
    dt = ScalarDT2D(maxeig,Problem.CFL,Mesh);
else
    dt = Problem.fixed_dt;
end

if(Output.show_plot)
    % set up rendering
    [TRI,xout,yout,interp] = GenInterpolators2D(Mesh.N, Mesh.N, Mesh.x, Mesh.y, Mesh.invV);
    Nq = 400;
    xmax = max(max(Mesh.x)); xmin = min(min(Mesh.x));
    ymax = max(max(Mesh.y)); ymin = min(min(Mesh.y));
    hx = (xmax-xmin)/Nq; hy = (ymax-ymin)/Nq;
    [xq,yq]=meshgrid(xmin:hx:xmax,ymin:hy:ymax);
end

% Initializing save arrays
Q_save     = cell(1,Problem.tstamps+2); % Save at approximated final time and actual final time
ind_save   = cell(1,Problem.tstamps+2); % Save at approximated final time and actual final time
ptc_hist   = [];
t_hist     = [];
Save_times = [linspace(0,Problem.FinalTime,Problem.tstamps+1),Problem.FinalTime];


tstep = 1; time = 0; stime = 1;

QG   = ApplyBCScalar2D(Q,time,Mesh);

% limit initial condition
ind  = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind,Limit.Limiter,Mesh);


perc_tcells = length(ind)*100/Mesh.K;
maxp_tcells = perc_tcells;
if(Output.show_plot)
    %plot solutions
    figure(1)
    plot(time,perc_tcells,'r.')
    xlim([0,Problem.FinalTime])
    ylim([0,100])
    xlabel('time')
    ylabel('Percentage of cells flagged')
    title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
    hold all
    drawnow
    
    figure(2);
    plot(xavg(ind), yavg(ind), 'r.');
    xlim(Output.xran)
    ylim(Output.yran)
    title(['Troubled cells t=',num2str(time)]);
    drawnow;
    
    figure(3)
    QC = Q(:,:,1);
    F  = scatteredInterpolant(Mesh.x(:),Mesh.y(:),QC(:)); interpc=F(xq,yq);
    contour(xq,yq,interpc,Output.clines,'k-','LineWidth',1);
    xlim(Output.xran)
    ylim(Output.yran)
    xlabel('x')
    ylabel('y')
    %PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
    
    figure(4)
    PlotField2D(Q(:,:,1),interp,TRI,xout,yout); xlim(Output.xran)
    ylim(Output.yran)
    drawnow;
    
    pause(0.5);
end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    ptc_hist(tstep)   = perc_tcells;
    t_hist(tstep)     = time;
    stime             = stime + 1;
end



% outer time step loop
while (time<Problem.FinalTime)
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % SSP-RK3 Stage 1
    rhsQ = ScalarRHS2D(Q, time,Problem,Mesh);
    Q1   = Q + dt*rhsQ;
    QG   = ApplyBCScalar2D(Q1,time,Mesh);
    ind1 = Find_Tcells2D(Q1(:,:,1),QG(:,:,1),Limit,Mesh,Net);
    Q1   = SlopeLimiter2D(Q1(:,:,1),QG(:,:,1),ind1,Limit.Limiter,Mesh);
    
    % SSP-RK3 Stage 2
    rhsQ = ScalarRHS2D(Q1, time,Problem,Mesh);
    Q2   = (3*Q + Q1 + dt*rhsQ)/4;
    QG   = ApplyBCScalar2D(Q2,time,Mesh);
    ind2 = Find_Tcells2D(Q2(:,:,1),QG(:,:,1),Limit,Mesh,Net);
    Q2   = SlopeLimiter2D(Q2(:,:,1),QG(:,:,1),ind2,Limit.Limiter,Mesh);
    
    % SSP-RK3 Stage 3
    rhsQ = ScalarRHS2D(Q2, time,Problem,Mesh);
    Q    = (Q + 2*Q2 + 2*dt*rhsQ)/3;
    QG   = ApplyBCScalar2D(Q,time,Mesh);
    ind3 = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
    Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind3,Limit.Limiter,Mesh);
    
    ind = unique([ind1,ind2, ind3]);
    
    % Increment time and compute new timestep
    time   = time+dt;
    
    
    
    
    perc_tcells = length(ind)*100/Mesh.K;
    maxp_tcells = max(maxp_tcells,perc_tcells);
    if(Output.show_plot)
        figure(1)
        plot(time,perc_tcells,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,100])
        xlabel('time')
        ylabel('Percentage of cells flagged')
        title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
        hold all
    end
    
    % Render every 25 time steps
    if(~mod(tstep,Output.plot_iter) || time == Problem.FinalTime)
        
        fprintf('   ---> Time = %f, dt = %f, number of time steps = %f\n',time,dt,tstep);
        
        if(Output.show_plot)
            figure(2);
            plot(xavg(ind), yavg(ind), 'r.');
            xlim(Output.xran)
            ylim(Output.yran)
            title(['Troubled cells t=',num2str(time)]);
            
            figure(3)
            QC = Q(:,:,1);
            F  = scatteredInterpolant(Mesh.x(:),Mesh.y(:),QC(:)); interpc=F(xq,yq);
            contour(xq,yq,interpc,Output.clines,'k-','LineWidth',1);
            xlim(Output.xran)
            ylim(Output.yran)
            xlabel('x')
            ylabel('y')
            
            
            figure(4)
            PlotField2D(Q(:,:,1),interp,TRI,xout,yout);
            xlim(Output.xran)
            ylim(Output.yran)
            drawnow;
            
            
            pause(0.1);
        end
        
    end
    
    if(Output.save_soln)
        ptc_hist(tstep+1) = perc_tcells;
        t_hist(tstep+1)   = time;
        if(abs(time-Save_times(stime))< 0.51*dt && stime <= length(Save_times)-1)
            Q_save{1,stime}   = Q;
            ind_save{1,stime} = ind;
            Save_times(stime) = time;
            stime             = stime + 1;
        end
    end
    
    maxeig = Get_scalar_eig2D(Q,Problem);
    if(isempty(Problem.fixed_dt))
        dt = ScalarDT2D(maxeig,Problem.CFL,Mesh);
    else
        dt = Problem.fixed_dt;
    end
    
    tstep = tstep+1;
end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
end

return;
