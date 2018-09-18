function [Q_save,ind_save,ptc_hist,t_hist,Save_times] = Scalar2D(Q,AdvectionVelocity,Save_soln)

% function Q = Scalar2D(Q)
% Purpose  : Integrate 2D Scalar equation using a SSP-RK3

Globals2D_DG;

% Get centroid
xavg = AVG2D*x;
yavg = AVG2D*y;

% Get max wave speed
maxeig = Get_scalar_eig2D(Q,model,AdvectionVelocity);

% compute initial timestep
if(isempty(fixed_dt))
    dt = ScalarDT2D(maxeig); 
else
    dt = fixed_dt;
end


Q_save     = cell(1,tstamps+2); % Save at approximated final time and actual final time
ind_save   = cell(1,tstamps+2); % Save at approximated final time and actual final time
ptc_hist   = [];
t_hist     = [];
Save_times = [linspace(0,FinalTime,tstamps+1),FinalTime];


tstep = 1; time = 0; stime = 1;

QG   = ApplyBCScalar2D(Q,time);

% limit initial condition
ind  = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter);
Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind,Limiter);

% set up rendering
[TRI,xout,yout,interp] = GenInterpolators2D(N, N, x, y, invV);

%plot solutions
figure(1)
perc_tcells = length(ind)*100/K;
maxp_tcells = perc_tcells;
plot(time,perc_tcells,'r.')
xlim([0,FinalTime])
ylim([0,100])
xlabel('time')
ylabel('Percentage of cells flagged')
title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(K),' cells'])
hold all
drawnow

figure(2);
plot(xavg(ind), yavg(ind), 'r.');
xlim([min(min(x)),max(max(x))])
ylim([min(min(y)),max(max(y))])
title(['Troubled cells t=',num2str(time)]);
drawnow;

figure(3)
PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;

% figure(4)
% trisurf(triplot, xplot, yplot, interpplot*Q(:,:,1));
% shading interp; view(2); title(['Solution, t=',num2str(time)]); axis tight; drawnow;

pause(0.5);

if(Save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    ptc_hist(tstep)   = perc_tcells;
    t_hist(tstep)     = time;
    stime             = stime + 1;
end



% outer time step loop
while (time<FinalTime)
            
    if(time+dt>FinalTime)
        dt = FinalTime-time;
    end
    
    % SSP-RK3 Stage 1
    rhsQ = ScalarRHS2D(Q, time,  AdvectionVelocity);
    Q1   = Q + dt*rhsQ;
    QG   = ApplyBCScalar2D(Q1,time);
    ind1 = Find_Tcells2D(Q1(:,:,1),QG(:,:,1),Indicator,Limiter);
    Q1   = SlopeLimiter2D(Q1(:,:,1),QG(:,:,1),ind1,Limiter);
    
    % SSP-RK3 Stage 2
    rhsQ = ScalarRHS2D(Q1, time,  AdvectionVelocity);
    Q2   = (3*Q + Q1 + dt*rhsQ)/4;
    QG   = ApplyBCScalar2D(Q2,time);
    ind2 = Find_Tcells2D(Q2(:,:,1),QG(:,:,1),Indicator,Limiter);
    Q2   = SlopeLimiter2D(Q2(:,:,1),QG(:,:,1),ind2,Limiter);
    
    % SSP-RK3 Stage 3
    rhsQ = ScalarRHS2D(Q2, time,  AdvectionVelocity);
    Q    = (Q + 2*Q2 + 2*dt*rhsQ)/3;
    QG   = ApplyBCScalar2D(Q,time);
    ind3 = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter);
    Q    = SlopeLimiter2D(Q(:,:,1),QG(:,:,1),ind3,Limiter);
    
    ind = unique([ind1,ind2, ind3]);
    
    % Increment time and compute new timestep
    time   = time+dt;
    
    
    
    figure(1)
    perc_tcells = length(ind)*100/K;
    
    maxp_tcells = max(maxp_tcells,perc_tcells);
    plot(time,perc_tcells,'r.')
    %xlim([0,FinalTime])
    ylim([0,100])
    xlabel('time')
    ylabel('Percentage of cells flagged')
    title(['Max = ',num2str(maxp_tcells),'%',' of ',num2str(K),' cells'])
    hold all
    
    % Render every 25 time steps
    if(~mod(tstep,100) || time == FinalTime)
        
        fprintf('   ---> Time = %f, dt = %f, number of time steps = %f\n',time,dt,tstep);
        
        figure(2);
        plot(xavg(ind), yavg(ind), 'r.');
        xlim([min(min(x)),max(max(x))])
        ylim([min(min(y)),max(max(y))])
        title(['Troubled cells t=',num2str(time)]);
        
        figure(3)
        PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
        
%         figure(4)
%         trisurf(triplot, xplot, yplot, interpplot*Q(:,:,1));
%         shading interp; view(2); title(['Solution, t=',num2str(time)]); axis tight;
        
        
        pause(0.1);
        
    end
    
    if(Save_soln)
        ptc_hist(tstep+1) = perc_tcells;
        t_hist(tstep+1)   = time;
        if(abs(time-Save_times(stime))< 0.51*dt && stime <= length(Save_times)-1)
            Q_save{1,stime}   = Q;
            ind_save{1,stime} = ind;
            Save_times(stime) = time;
            stime             = stime + 1;
        end
    end
    
    maxeig = Get_scalar_eig2D(Q,model,AdvectionVelocity);
    if(isempty(fixed_dt))
        dt = ScalarDT2D(maxeig);
    else
        dt = fixed_dt;
    end
    
    tstep = tstep+1;
end

if(Save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
end

return;
