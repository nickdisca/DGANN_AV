function [Q_save,ind_save,visc_save,ptc_hist,pnc_hist,maxvisc_hist,t_hist,Save_times] = Euler2D(Q,Problem,Mesh,Limit,Net,Viscosity,NetVisc,Output)

% Get centroid
xavg = Mesh.AVG2D*Mesh.x;
yavg = Mesh.AVG2D*Mesh.y;


% % compute initial timestep
% if(isempty(Problem.fixed_dt))
%     dt      = EulerDT2D(Q,Problem,Mesh);
%     CFL     = Problem.CFL;
% else
%     dt      = Problem.fixed_dt;
%     CFL     = EulerCFL2D(Q,dt,Problem,Mesh);
% end

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
pnc_hist   = [];
maxvisc_hist  = []; 
t_hist     = [];
Save_times = [linspace(0,Problem.FinalTime,Problem.tstamps+1),Problem.FinalTime];

dt = 0; tstep = 0; time = 0; stime = 1;



% limit initial condition
QG          = ApplyBCEuler2D(Q,time,Problem,Mesh);
ind         = Tcells_Euler_type2D(Q,QG,Problem.gas_gamma,Limit,Mesh,Net);
Q           = SlopeLimit_Euler_type2D(Q,QG,ind,Problem.gas_gamma,Limit,Mesh);
[Q,ind_neg] = Euler_positivity_fix2D(Q,Problem.gas_gamma,Mesh.AVG2D);


perc_tcells  = length(ind)*100/Mesh.K;
maxp_tcells  = perc_tcells;
perc_tcells2 = length(ind_neg)*100/Mesh.K;
maxp_tcells2 = perc_tcells2;
if(Output.show_plot)
    figure(1); subplot(3,1,1);
    plot(time,perc_tcells,'r.')
    xlim([0,Problem.FinalTime])
    ylim([0,100])
    xlabel({'$t$'},'interpreter','latex')
    ylabel('Percentage of cells flagged')
    title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
    hold all
    drawnow
    
    figure(1); subplot(3,1,2);
    plot(time,perc_tcells2,'r.')
    xlim([0,Problem.FinalTime])
    ylim([0,100])
    xlabel('time')
    ylabel('Percentage of cells fixed due to loss of positivity')
    title(['Max = ',num2str(maxp_tcells2),'%',' of ',num2str(Mesh.K),' cells'])
    hold all
    drawnow
    
    
    figure(3);
    plot(xavg(ind), yavg(ind), 'r.');
    hold all
    plot(xavg(ind_neg), yavg(ind_neg), 'bx');
    xlim(Output.xran)
    ylim(Output.yran)
    title(['Troubled cells (red) and neg-fixed cells (blue) t=',num2str(time)]);
    drawnow;
    hold off
    
    PlotEulerFields2D(Q,time,Problem.gas_gamma,Problem.gas_const,Mesh.x,Mesh.y,xq,yq,TRI,xout,yout,interp,Output,2,100);
    
    pause(0.1);
end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    ptc_hist(tstep+1)   = perc_tcells;
    pnc_hist(tstep+1)   = perc_tcells2;
    t_hist(tstep+1)     = time;
    stime             = stime + 1;
end

% Temporary storage for solution at previous tstep (needed for EV)
Q_tmp=NaN*zeros(size(Q,1)*size(Q,2)*size(Q,3),2); Q_tmp(:,2)=Q(:);


% outer time step loop
while (time<Problem.FinalTime)

    % Get max wave speed
    vel_u = Q(:,:,2)./Q(:,:,1); vel_v = Q(:,:,3)./Q(:,:,1);
    pres = (Problem.gas_gamma-1.0)*(Q(:,:,4) - Q(:,:,1).*(vel_u.^2+vel_v.^2)/2); 
    c_sound = sqrt(abs(Problem.gas_gamma*pres./Q(:,:,1)));
    lambda = (sqrt (vel_u.^2 + vel_v.^2) + c_sound);
    
    %Compute artificial viscosity
    Qold = reshape(Q_tmp(:,1),Mesh.Np,Mesh.K,4);
    mu_piece = Euler2D_viscosity(Q, Qold, lambda, Viscosity, Problem, Mesh, dt, tstep, NetVisc);
    mu_piece=max(mu_piece,0);
    mu_vals=Scalar2D_smooth_viscosity(mu_piece,Mesh); 
    mu_vals=max(mu_vals,0);
    maxvisc=max(abs(mu_vals(:)));
    
    %Apply same AV for all equations
    mu_vals=repmat(mu_vals,1,1,4);
    
    if(isempty(Problem.fixed_dt))
        dt = EulerDT2D(lambda,maxvisc,Problem.CFL,Mesh);
    else
        dt = Problem.fixed_dt;
    end
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    if strcmp(Problem.RK,'SSP3')
        
        % SSP-RK3 Stage 1
        rhsQ = EulerRHS2D(Q, time, mu_vals, Problem,Mesh);
        Q1   = Q + dt*rhsQ;
        
        QG              = ApplyBCEuler2D(Q1,time,Problem,Mesh);
        ind1            = Tcells_Euler_type2D(Q1,QG,Problem.gas_gamma,Limit,Mesh,Net);
        Q1              = SlopeLimit_Euler_type2D(Q1,QG,ind1,Problem.gas_gamma,Limit,Mesh);
        [Q1,ind_neg1]   = Euler_positivity_fix2D(Q1,Problem.gas_gamma,Mesh.AVG2D);
        %assert(isempty(find(Q1(:,:,1) <= 0,1)));
        %assert(isempty(find(Euler_Pressure2D(Q1,Problem.gas_gamma) <= 0,1)));
        
        
        % SSP-RK3 Stage 2
        rhsQ = EulerRHS2D(Q1, time + dt, mu_vals, Problem,Mesh);
        Q2   = (3*Q + Q1 + dt*rhsQ)/4;
        
        QG              = ApplyBCEuler2D(Q2,time+dt,Problem,Mesh);
        ind2            = Tcells_Euler_type2D(Q2,QG,Problem.gas_gamma,Limit,Mesh,Net);
        Q2              = SlopeLimit_Euler_type2D(Q2,QG,ind2,Problem.gas_gamma,Limit,Mesh);
        [Q2,ind_neg2]   = Euler_positivity_fix2D(Q2,Problem.gas_gamma,Mesh.AVG2D);
        %assert(isempty(find(Q2(:,:,1) <= 0,1)));
        %assert(isempty(find(Euler_Pressure2D(Q2,Problem.gas_gamma) <= 0,1)));
        
        % SSP-RK3 Stage 3
        rhsQ = EulerRHS2D(Q2, time + 0.5*dt, mu_vals, Problem,Mesh);
        Q    = (Q + 2*Q2 + 2*dt*rhsQ)/3;
        
        QG             = ApplyBCEuler2D(Q,time+ 0.5*dt,Problem,Mesh);
        ind3           = Tcells_Euler_type2D(Q,QG,Problem.gas_gamma,Limit,Mesh,Net);
        Q              = SlopeLimit_Euler_type2D(Q,QG,ind3,Problem.gas_gamma,Limit,Mesh);
        [Q,ind_neg3]   = Euler_positivity_fix2D(Q,Problem.gas_gamma,Mesh.AVG2D);
        %assert(isempty(find(Q(:,:,1) <= 0,1)));
        %assert(isempty(find(Euler_Pressure2D(Q,Problem.gas_gamma) <= 0,1)));
        
        ind = unique([ind1,ind2, ind3]);
        ind_neg = unique([ind_neg1,ind_neg2, ind_neg3]);
        
    elseif strcmp(Problem.RK,'LS54')
        
        A_RK(1)=0; A_RK(2)=-0.4178904745; A_RK(3)=-1.192151694643; A_RK(4)=-1.697784692471; A_RK(5)=-1.514183444257;
        B_RK(1)=0.1496590219993; B_RK(2)=0.3792103129999; B_RK(3)=0.8229550293869; B_RK(4)=0.6994504559488; B_RK(5)=0.1530572479681;
        c_RK(1)=0; c_RK(2)=0.1496590219993; c_RK(3)=0.3704009573644; c_RK(4)=0.6222557631345; c_RK(5)=0.9582821306748;
        nsteps=5;
        
        VV=zeros(size(Q));
        ind=[];
        ind_neg=[];
        
        for index=1:nsteps
            
            rhsQ = EulerRHS2D(Q, time + c_RK(index)*dt, mu_vals, Problem,Mesh);
            
            VV=A_RK(index)*VV+dt*rhsQ;
            
            Q=Q+B_RK(index)*VV;
            
            QG   = ApplyBCEuler2D(Q,time+ c_RK(index)*dt,Problem,Mesh);
            ind_substage  = Tcells_Euler_type2D(Q,QG,Problem.gas_gamma,Limit,Mesh,Net);
            Q = SlopeLimit_Euler_type2D(Q,QG,ind_substage,Problem.gas_gamma,Limit,Mesh);
            [Q,ind_neg_substage] = Euler_positivity_fix2D(Q,Problem.gas_gamma,Mesh.AVG2D);
            
            ind = unique([ind,ind_substage]);
            ind_neg = unique([ind_neg,ind_neg_substage]);
        end
        
    else
        
        error('Time integration scheme not defined');
    end
    
    
    % Increment time and timestep
    time   = time+dt;
    tstep  = tstep +1;
    
    % Increment saved variables (for EV)
    Q_tmp(:,1)=Q_tmp(:,2); Q_tmp(:,2)=Q(:);
    
    perc_tcells = length(ind)*100/Mesh.K;
    maxp_tcells = max(maxp_tcells,perc_tcells);
    perc_tcells2 = length(ind_neg)*100/Mesh.K;
    maxp_tcells2 = max(maxp_tcells2,perc_tcells2);
    if(Output.show_plot)
        
        figure(1); subplot(3,1,1);
        plot(time,perc_tcells,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,100])
        xlabel('time')
        ylabel('Percentage of cells flagged')
        title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
        hold all
        
        figure(1); subplot(3,1,2);
        plot(time,perc_tcells2,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,100])
        xlabel('time')
        ylabel('Percentage of cells fixed due to loss of positivity')
        title(['Max = ',num2str(maxp_tcells2),'%',' of ',num2str(Mesh.K),' cells'])
        hold all
        
        figure(1); subplot(3,1,3);
        plot(time,maxvisc,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,0.2])
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
    if(~mod(tstep,Output.plot_iter) || time == Problem.FinalTime)
        
        fprintf('   ---> Time = %f, CFL = %f, dt = %f, number of time steps = %d\n',time,Problem.CFL,dt,tstep);
        
        if(Output.show_plot)
            figure(3);
            plot(xavg(ind), yavg(ind), 'r.');
            hold all
            plot(xavg(ind_neg), yavg(ind_neg), 'bx');
            xlim(Output.xran)
            ylim(Output.yran)
            title(['Troubled cells (red) and neg-fixed cells (blue) t=',num2str(time)]);
            drawnow;
            hold off
            
            PlotEulerFields2D(Q,time,Problem.gas_gamma,Problem.gas_const,Mesh.x,Mesh.y,xq,yq,TRI,xout,yout,interp,Output,2,100);
            
            figure(4);
            mu_plot=mu_vals(:,:,1);
            F = scatteredInterpolant(Mesh.x(:),Mesh.y(:),mu_plot(:)); interpc=F(xq,yq);
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
        pnc_hist(tstep+1) = perc_tcells2;
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
