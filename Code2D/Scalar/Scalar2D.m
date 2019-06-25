function [Q_save,ind_save,visc_save,ptc_hist,maxvisc_hist,t_hist,Save_times] = Scalar2D(Q,Problem,Mesh,Limit,Net,Viscosity,NetVisc,Output)

% function Q = Scalar2D(Q)
% Purpose  : Integrate 2D Scalar equation using a RK scheme


% Get centroid
xavg = Mesh.AVG2D*Mesh.x;
yavg = Mesh.AVG2D*Mesh.y;

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
visc_save   = cell(1,Problem.tstamps+2); % Save at approximated final time and actual final time
ptc_hist   = [];
maxvisc_hist  = []; 
t_hist     = [];
Save_times = [linspace(0,Problem.FinalTime,Problem.tstamps+1),Problem.FinalTime];


dt = 0; tstep = 0; time = 0; stime = 1;

QG   = ApplyBCScalar2D(Q,time,Mesh);

% limit initial condition
ind  = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind,Limit.Limiter,Mesh);


perc_tcells = length(ind)*100/Mesh.K;
maxp_tcells = perc_tcells;
if(Output.show_plot)
    %plot solutions
    figure(1); subplot(2,1,1);
    plot(time,perc_tcells,'r.')
    xlim([0,Problem.FinalTime])
    ylim([0,100])
    xlabel({'$t$'},'interpreter','latex')
    ylabel('Percentage of cells flagged')
    title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
    hold all
    drawnow
    
    figure(3);
    plot(xavg(ind), yavg(ind), 'r.');
    xlim(Output.xran)
    ylim(Output.yran)
    title(['Troubled cells t=',num2str(time)]);
    drawnow;
    
    figure(2);
    QC = Q(:,:,1);
    F  = scatteredInterpolant(Mesh.x(:),Mesh.y(:),QC(:)); interpc=F(xq,yq);
    contourf(xq,yq,interpc,Output.clines,'k-','LineWidth',1);
    colormap jet;
    xlim(Output.xran)
    ylim(Output.yran)
    xlabel({'$x$'},'interpreter','latex')
    ylabel({'$y$'},'interpreter','latex')
    %PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
    
    pause(0.1);
end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    ptc_hist(tstep+1)   = perc_tcells;
    t_hist(tstep+1)     = time;
    stime             = stime + 1;
end

% Temporary storage for solution at previous tstep (needed for EV)
Q_tmp=NaN*zeros(size(Q,1)*size(Q,2)*size(Q,3),2); Q_tmp(:,2)=Q(:);



% outer time step loop
while (time<Problem.FinalTime)
    
    % Get max wave speed
    maxeig = Get_scalar_eig2D(Q,Problem);
    
    % Compute artificial viscosity
    Qold = reshape(Q_tmp(:,1),Mesh.Np,Mesh.K);
    mu_piece = Scalar2D_viscosity(Q, Qold, maxeig, Viscosity, Problem, Mesh, dt, tstep, NetVisc);
    mu_piece=max(mu_piece,0);
    mu_vals=Scalar2D_smooth_viscosity(mu_piece,Mesh); 
    mu_vals=max(mu_vals,0);
    maxvisc=max(abs(mu_vals(:)));
    
    % compute timestep
    if(isempty(Problem.fixed_dt))
        dt = ScalarDT2D(maxeig,maxvisc,Problem.CFL,Mesh);
    else
        dt = Problem.fixed_dt;
    end
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % Save aux variable (needed for EV)
    Q_tmp(:,2)=Q(:);
    
    % 3rd order SSP Runge-Kutta
    if strcmp(Problem.RK,'SSP3')
    
        % SSP-RK3 Stage 1
        rhsQ = ScalarRHS2D(Q, time,mu_vals,Problem,Mesh);
        Q1   = Q + dt*rhsQ;
        QG   = ApplyBCScalar2D(Q1,time,Mesh);
        ind1 = Find_Tcells2D(Q1(:,:,1),QG(:,:,1),Limit,Mesh,Net);
        Q1   = SlopeLimiter2D(Q1(:,:,1),QG(:,:,1),ind1,Limit.Limiter,Mesh);
        
        % SSP-RK3 Stage 2
        rhsQ = ScalarRHS2D(Q1, time,mu_vals,Problem,Mesh);
        Q2   = (3*Q + Q1 + dt*rhsQ)/4;
        QG   = ApplyBCScalar2D(Q2,time,Mesh);
        ind2 = Find_Tcells2D(Q2(:,:,1),QG(:,:,1),Limit,Mesh,Net);
        Q2   = SlopeLimiter2D(Q2(:,:,1),QG(:,:,1),ind2,Limit.Limiter,Mesh);
        
        % SSP-RK3 Stage 3
        rhsQ = ScalarRHS2D(Q2, time,mu_vals,Problem,Mesh);
        Q    = (Q + 2*Q2 + 2*dt*rhsQ)/3;
        QG   = ApplyBCScalar2D(Q,time,Mesh);
        ind3 = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
        Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind3,Limit.Limiter,Mesh);
        
        ind = unique([ind1,ind2, ind3]);
    
    elseif strcmp(Problem.RK,'LS54')
        
        A_RK(1)=0; A_RK(2)=-0.4178904745; A_RK(3)=-1.192151694643; A_RK(4)=-1.697784692471; A_RK(5)=-1.514183444257;
        B_RK(1)=0.1496590219993; B_RK(2)=0.3792103129999; B_RK(3)=0.8229550293869; B_RK(4)=0.6994504559488; B_RK(5)=0.1530572479681;
        c_RK(1)=0; c_RK(2)=0.1496590219993; c_RK(3)=0.3704009573644; c_RK(4)=0.6222557631345; c_RK(5)=0.9582821306748;
        nsteps=5;
        
        VV=zeros(size(Q));
        ind=[];
        
        for index=1:nsteps
            
            rhsQ = ScalarRHS2D(Q, time,mu_vals,Problem,Mesh);
            
            VV=A_RK(index)*VV+dt*rhsQ;
            
            Q=Q+B_RK(index)*VV;
            
            QG   = ApplyBCScalar2D(Q,time,Mesh);
            ind_substage = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
            Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind_substage,Limit.Limiter,Mesh);
        
            ind = unique([ind,ind_substage]);
        end
        
    else
        
        error('Time integration scheme not defined');
    end
    
    % Increment time and timestep
    time   = time+dt;
    tstep  = tstep +1;
    
    % Increment saved variables
    Q_tmp(:,1)=Q_tmp(:,2); Q_tmp(:,2)=Q(:);
    
    perc_tcells = length(ind)*100/Mesh.K;
    maxp_tcells = max(maxp_tcells,perc_tcells);
    if(Output.show_plot)
        figure(1); subplot(2,1,1);
        plot(time,perc_tcells,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,100])
        xlabel('time')
        ylabel('Percentage of cells flagged')
        title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
        hold all
        
        figure(1); subplot(2,1,2);
        plot(time,maxvisc,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,0.1])
        xlabel('time')
        ylabel('Maximum value of artificial viscosity')
        title('Max AV')
        hold all
    end
    
    % Check that solution is not diverging
    max_allowed_val=100;
    if max(abs(Q(:)))>=max_allowed_val
        error('Solution is diverging!');
    end
    
    % Render every plot_iter time steps
    if(~mod(tstep+1,Output.plot_iter) || time == Problem.FinalTime)
        
        fprintf('   ---> Time = %f, dt = %f, number of time steps = %d\n',time,dt,tstep+1);
        
        if(Output.show_plot)
            figure(3);
            plot(xavg(ind), yavg(ind), 'r.');
            xlim(Output.xran)
            ylim(Output.yran)
            title(['Troubled cells t=',num2str(time)]);
            
            figure(2)
            QC = Q(:,:,1);
            F  = scatteredInterpolant(Mesh.x(:),Mesh.y(:),QC(:)); interpc=F(xq,yq);
            contourf(xq,yq,interpc,Output.clines,'LineWidth',1);
            colormap jet;
            xlim(Output.xran)
            ylim(Output.yran)
            xlabel({'$x$'},'interpreter','latex')
            ylabel({'$y$'},'interpreter','latex')
            
            figure(4)
            F  = scatteredInterpolant(Mesh.x(:),Mesh.y(:),mu_vals(:)); interpc=F(xq,yq);
            contourf(xq,yq,interpc,'LineWidth',1);
            colormap jet;
            xlabel({'$x$'},'interpreter','latex')
            ylabel({'$y$'},'interpreter','latex') 
            title(['Viscosity t=',num2str(time)]);
            
            
            pause(0.1);
        end
        
    end
    
    if(Output.save_soln)
        ptc_hist(tstep+1) = perc_tcells;
        maxvisc_hist(tstep+1) = maxvisc;
        t_hist(tstep+1)   = time;
        if(abs(time-Save_times(stime))< 0.51*dt && stime <= length(Save_times)-1)
            Q_save{1,stime}   = Q;
            ind_save{1,stime} = ind;
            visc_save{1,stime} = mu_vals;
            Save_times(stime) = time;
            stime             = stime + 1;
        end
    end
    

end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    visc_save{1,stime} = mu_vals;
end

return;
