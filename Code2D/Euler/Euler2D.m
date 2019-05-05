function [Q_save,ind_save,ptc_hist,pnc_hist,t_hist,Save_times] = Euler2D(Q,Problem,Mesh,Limit,Net,Output)

% Get centroid
xavg = Mesh.AVG2D*Mesh.x;
yavg = Mesh.AVG2D*Mesh.y;


% compute initial timestep
if(isempty(Problem.fixed_dt))
    dt      = EulerDT2D(Q,Problem,Mesh);
    CFL     = Problem.CFL;
else
    dt      = Problem.fixed_dt;
    CFL     = EulerCFL2D(Q,dt,Problem,Mesh);
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
pnc_hist   = [];
t_hist     = [];
Save_times = [linspace(0,Problem.FinalTime,Problem.tstamps+1),Problem.FinalTime];

tstep = 1; time = 0; stime = 1;

% set up rendering
% [TRI,xout,yout,interp] = GenInterpolators2D(N, N, x, y, invV);
% [TRIG,xGout,yGout,interpG] = GenInterpolators2D(N, N, xG, yG, invV);

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
    figure(1)
    plot(time,perc_tcells,'r.')
    xlim([0,Problem.FinalTime])
    ylim([0,100])
    xlabel('time')
    ylabel('Percentage of cells flagged')
    title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
    hold all
    drawnow
    
    figure(2)
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
    
    PlotEulerFields2D(Q,time,Problem.gas_gamma,Problem.gas_const,Mesh.x,Mesh.y,xq,yq,TRI,xout,yout,interp,Output,4,5);
    
    pause(0.1);
    
end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    ptc_hist(tstep)   = perc_tcells;
    pnc_hist(tstep)   = perc_tcells2;
    t_hist(tstep)     = time;
    stime             = stime + 1;
end


%resQ = zeros(Np,K);

% outer time step loop
while (time<Problem.FinalTime)
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % SSP-RK3 Stage 1
    rhsQ = EulerRHS2D(Q, time, Problem,Mesh);
    Q1   = Q + dt*rhsQ;
    
    QG              = ApplyBCEuler2D(Q1,time,Problem,Mesh);
    ind1            = Tcells_Euler_type2D(Q1,QG,Problem.gas_gamma,Limit,Mesh,Net);
    Q1              = SlopeLimit_Euler_type2D(Q1,QG,ind1,Problem.gas_gamma,Limit,Mesh);
    [Q1,ind_neg1]   = Euler_positivity_fix2D(Q1,Problem.gas_gamma,Mesh.AVG2D);
    assert(isempty(find(Q1(:,:,1) <= 0,1)));
    assert(isempty(find(Euler_Pressure2D(Q1,Problem.gas_gamma) <= 0,1)));
    
    
    % SSP-RK3 Stage 2
    rhsQ = EulerRHS2D(Q1, time + dt, Problem,Mesh);
    Q2   = (3*Q + Q1 + dt*rhsQ)/4;
    
    QG              = ApplyBCEuler2D(Q2,time+dt,Problem,Mesh);
    ind2            = Tcells_Euler_type2D(Q2,QG,Problem.gas_gamma,Limit,Mesh,Net);
    Q2              = SlopeLimit_Euler_type2D(Q2,QG,ind2,Problem.gas_gamma,Limit,Mesh);
    [Q2,ind_neg2]   = Euler_positivity_fix2D(Q2,Problem.gas_gamma,Mesh.AVG2D);
    assert(isempty(find(Q2(:,:,1) <= 0,1)));
    assert(isempty(find(Euler_Pressure2D(Q2,Problem.gas_gamma) <= 0,1)));
    
    % SSP-RK3 Stage 3
    rhsQ = EulerRHS2D(Q2, time + 0.5*dt, Problem,Mesh);
    Q    = (Q + 2*Q2 + 2*dt*rhsQ)/3;
    
    QG             = ApplyBCEuler2D(Q,time+ 0.5*dt,Problem,Mesh);
    ind3           = Tcells_Euler_type2D(Q,QG,Problem.gas_gamma,Limit,Mesh,Net);
    Q              = SlopeLimit_Euler_type2D(Q,QG,ind3,Problem.gas_gamma,Limit,Mesh);
    [Q,ind_neg3]   = Euler_positivity_fix2D(Q,Problem.gas_gamma,Mesh.AVG2D);
    assert(isempty(find(Q(:,:,1) <= 0,1)));
    assert(isempty(find(Euler_Pressure2D(Q,Problem.gas_gamma) <= 0,1)));
    
    ind = unique([ind1,ind2, ind3]);
    ind_neg = unique([ind_neg1,ind_neg2, ind_neg3]);
    
    %     ind = [];
    %     for INTRK = 1:5
    %         rhsQ  = EulerRHS2D(Q, time, gas_gamma, gas_const);
    %         resQ  = rk4a(INTRK)*resQ + dt*rhsQ;
    %         Q     = Q + rk4b(INTRK)*resQ;
    %         indrk = Tcells_Euler_type2D(Q,gas_gamma, gas_const,Indicator,ind_var,Limiter);
    %         Q     = SlopeLimit_Euler_type2D(Q,gas_gamma,gas_const,Limiter,lim_var,indrk);
    %         assert(isempty(find(Q(:,:,1) <= 0)));
    %         assert(isempty(find(Euler_Pressure2D(Q,gas_gamma) <= 0)));
    %         ind   = [ind,indrk];
    %     end
    %     ind = unique(ind);
    
    
    
    % Increment time and compute new timestep
    time   = time+dt;
    
    perc_tcells = length(ind)*100/Mesh.K;
    maxp_tcells = max(maxp_tcells,perc_tcells);
    perc_tcells2 = length(ind_neg)*100/Mesh.K;
    maxp_tcells2 = max(maxp_tcells2,perc_tcells2);
    if(Output.show_plot)
        figure(1)
        plot(time,perc_tcells,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,100])
        xlabel('time')
        ylabel('Percentage of cells flagged')
        title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(Mesh.K),' cells'])
        hold all
        
        figure(2)
        
        plot(time,perc_tcells2,'r.')
        xlim([0,Problem.FinalTime])
        ylim([0,100])
        xlabel('time')
        ylabel('Percentage of cells fixed due to loss of positivity')
        title(['Max = ',num2str(maxp_tcells2),'%',' of ',num2str(Mesh.K),' cells'])
        hold all
    end
    
    % Render every 25 time steps
    if(~mod(tstep,Output.plot_iter) || time == Problem.FinalTime)
        
        fprintf('   ---> Time = %f, CFL = %f, dt = %f, number of time steps = %f\n',time,CFL,dt,tstep);
        
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
            
            PlotEulerFields2D(Q,time,Problem.gas_gamma,Problem.gas_const,Mesh.x,Mesh.y,xq,yq,TRI,xout,yout,interp,Output,4,5);
            
            pause(0.1);
        end
        
    end
    
    if(Output.save_soln)
        ptc_hist(tstep+1) = perc_tcells;
        pnc_hist(tstep+1) = perc_tcells2;
        t_hist(tstep+1)   = time;
        if(abs(time-Save_times(stime))< 0.51*dt && stime <= length(Save_times)-1)
            Q_save{1,stime}   = Q;
            ind_save{1,stime} = ind;
            Save_times(stime) = time;
            stime             = stime + 1;
        end
    end
    
    if(isempty(Problem.fixed_dt))
        dt      = EulerDT2D(Q,Problem,Mesh);
    else
        dt      = Problem.fixed_dt;
        CFL     = EulerCFL2D(Q,dt,Problem,Mesh);
    end
    
    tstep = tstep+1;
end

if(Output.save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
end

return;
