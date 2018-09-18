function [Q_save,ind_save,ptc_hist,pnc_hist,t_hist,Save_times] = Euler2D(Q,gas_gamma,gas_const,Save_soln)

% function Q = Euler2D(Q)
% Purpose  : Integrate 2D Euler equation using a SSP-RK3

Globals2D_DG;

% Get centroid
xavg = AVG2D*x;
yavg = AVG2D*y;

% compute initial timestep
if(isempty(fixed_dt))
    dt = EulerDT2D(Q,gas_gamma);
else
    dt  = fixed_dt;
    CFL = EulerCFL2D(Q,dt,gas_gamma);
end


Q_save     = cell(1,tstamps+2); % Save at approximated final time and actual final time
ind_save   = cell(1,tstamps+2); % Save at approximated final time and actual final time
ptc_hist   = [];
pnc_hist   = [];
t_hist     = [];
Save_times = [linspace(0,FinalTime,tstamps+1),FinalTime];

tstep = 1; time = 0; stime = 1;

% set up rendering
[TRI,xout,yout,interp] = GenInterpolators2D(N, N, x, y, invV);
[TRIG,xGout,yGout,interpG] = GenInterpolators2D(N, N, xG, yG, invV);

% limit initial condition
QG          = ApplyBCEuler2D(Q,time,gas_gamma);
% figure(3)
% PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
% hold all
% figure(3)
% PlotField2D(QG(:,:,1),interpG,TRIG,xGout,yGout); axis tight; drawnow;

ind         = Tcells_Euler_type2D(Q,QG,gas_gamma, gas_const,Indicator,ind_var,Limiter);
Q           = SlopeLimit_Euler_type2D(Q,QG,gas_gamma,gas_const,Limiter,lim_var,ind);
[Q,ind_neg] = Euler_positivity_fix2D(Q,gas_gamma,AVG2D);




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
hold all
plot(xavg(ind_neg), yavg(ind_neg), 'bx');
xlim([min(min(x)),max(max(x))])
ylim([min(min(y)),max(max(y))])
title(['Troubled cells t=',num2str(time)]);
drawnow;
hold off

figure(3)
PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;

figure(4)
PlotField2D(Euler_Pressure2D(Q,gas_gamma),interp,TRI,xout,yout); axis tight; drawnow;


figure(5)
perc_tcells2 = length(ind_neg)*100/K;
maxp_tcells2 = perc_tcells2;
plot(time,perc_tcells2,'r.')
xlim([0,FinalTime])
ylim([0,100])
xlabel('time')
ylabel('Percentage of cells flagged')
title(['Max = ',num2str(maxp_tcells2),'%',' of ',num2str(K),' cells'])
hold all
drawnow

pause(0.5);

if(Save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
    ptc_hist(tstep)   = perc_tcells;
    pnc_hist(tstep)   = perc_tcells2;
    t_hist(tstep)     = time;
    stime             = stime + 1;
end


%resQ = zeros(Np,K);

% outer time step loop
while (time<FinalTime)
    
    if(time+dt>FinalTime)
        dt = FinalTime-time;
    end
    
    % SSP-RK3 Stage 1
    rhsQ = EulerRHS2D(Q, time, gas_gamma, gas_const);
    Q1   = Q + dt*rhsQ;
    %rhoprange = [min(min(Q1(:,:,1))), max(max(Q1(:,:,1))), min(min(Euler_Pressure2D(Q1,gas_gamma))), max(max(Euler_Pressure2D(Q1,gas_gamma)))]
    QG              = ApplyBCEuler2D(Q1,time,gas_gamma);
    ind1            = Tcells_Euler_type2D(Q1,QG,gas_gamma, gas_const,Indicator,ind_var,Limiter);
    Q1              = SlopeLimit_Euler_type2D(Q1,QG,gas_gamma,gas_const,Limiter,lim_var,ind1);
    [Q1,ind_neg1]   = Euler_positivity_fix2D(Q1,gas_gamma,AVG2D);
    %rhoprange = [min(min(Q1(:,:,1))), max(max(Q1(:,:,1))), min(min(Euler_Pressure2D(Q1,gas_gamma))), max(max(Euler_Pressure2D(Q1,gas_gamma)))]
    
%         figure(2);
%     plot(xavg(ind1), yavg(ind1), 'r.');
%     xlim([min(min(x)),max(max(x))])
%     ylim([min(min(y)),max(max(y))])
%     title(['Troubled cells t=',num2str(time)]);
%     drawnow;
%         figure(3)
%     PlotField2D(Q1(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
%     
%     figure(4)
%     PlotField2D(Euler_Pressure2D(Q1,gas_gamma),interp,TRI,xout,yout); axis tight; drawnow;
    
    assert(isempty(find(Q1(:,:,1) <= 0)));
    %     x(Euler_Pressure2D(Q1,gas_gamma) <= 0)
    %     y(Euler_Pressure2D(Q1,gas_gamma) <= 0)
    %     iii = find(Euler_Pressure2D(Q1,gas_gamma) <= 0)
    %     figure(2)
    %     hold all
    %     plot(x(iii), y(iii), 'bx');
    %     xlim([min(min(x)),max(max(x))])
    % ylim([min(min(y)),max(max(y))])
    %    PlotMesh2D();
    % hold off
    assert(isempty(find(Euler_Pressure2D(Q1,gas_gamma) <= 0)));
    
    
    % SSP-RK3 Stage 2
    rhsQ = EulerRHS2D(Q1, time + dt, gas_gamma, gas_const);
    Q2   = (3*Q + Q1 + dt*rhsQ)/4;
    %rhoprange = [min(min(Q2(:,:,1))), max(max(Q2(:,:,1))), min(min(Euler_Pressure2D(Q2,gas_gamma))), max(max(Euler_Pressure2D(Q2,gas_gamma)))]
    QG              = ApplyBCEuler2D(Q2,time+ dt,gas_gamma);
    ind2            = Tcells_Euler_type2D(Q2,QG,gas_gamma, gas_const,Indicator,ind_var,Limiter);
    Q2              = SlopeLimit_Euler_type2D(Q2,QG,gas_gamma,gas_const,Limiter,lim_var,ind2);
    [Q2,ind_neg2]   = Euler_positivity_fix2D(Q2,gas_gamma,AVG2D);
    %rhoprange = [min(min(Q2(:,:,1))), max(max(Q2(:,:,1))), min(min(Euler_Pressure2D(Q2,gas_gamma))), max(max(Euler_Pressure2D(Q2,gas_gamma)))]
    %     figure(2);
    % plot(xavg(ind2), yavg(ind2), 'r.');
    % xlim([min(min(x)),max(max(x))])
    % ylim([min(min(y)),max(max(y))])
    % title(['Troubled cells t=',num2str(time)]);
    % drawnow;
    %     figure(3)
    % PlotField2D(Q2(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
    %
    % figure(4)
    % PlotField2D(Euler_Pressure2D(Q2,gas_gamma),interp,TRI,xout,yout); axis tight; drawnow;
    assert(isempty(find(Q2(:,:,1) <= 0)));
    %     x(Euler_Pressure2D(Q2,gas_gamma) <= 0)
    %     y(Euler_Pressure2D(Q2,gas_gamma) <= 0)
    assert(isempty(find(Euler_Pressure2D(Q2,gas_gamma) <= 0)));
    
    % SSP-RK3 Stage 3
    rhsQ = EulerRHS2D(Q2, time + 0.5*dt, gas_gamma, gas_const);
    Q    = (Q + 2*Q2 + 2*dt*rhsQ)/3;
    %rhoprange = [min(min(Q(:,:,1))), max(max(Q(:,:,1))), min(min(Euler_Pressure2D(Q,gas_gamma))), max(max(Euler_Pressure2D(Q,gas_gamma)))]
    QG             = ApplyBCEuler2D(Q,time+ 0.5*dt,gas_gamma);
    ind3           = Tcells_Euler_type2D(Q,QG,gas_gamma, gas_const,Indicator,ind_var,Limiter);
    Q              = SlopeLimit_Euler_type2D(Q,QG,gas_gamma,gas_const,Limiter,lim_var,ind3);
    [Q,ind_neg3]   = Euler_positivity_fix2D(Q,gas_gamma,AVG2D);
    %rhoprange = [min(min(Q(:,:,1))), max(max(Q(:,:,1))), min(min(Euler_Pressure2D(Q,gas_gamma))), max(max(Euler_Pressure2D(Q,gas_gamma)))]
    
    %     figure(2);
    % plot(xavg(ind3), yavg(ind3), 'r.');
    % xlim([min(min(x)),max(max(x))])
    % ylim([min(min(y)),max(max(y))])
    % title(['Troubled cells t=',num2str(time)]);
    % drawnow;
    %     figure(3)
    % PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
    %
    % figure(4)
    % PlotField2D(Euler_Pressure2D(Q,gas_gamma),interp,TRI,xout,yout); axis tight; drawnow;
    assert(isempty(find(Q(:,:,1) <= 0)));
    %     x(Euler_Pressure2D(Q,gas_gamma) <= 0)
    %     y(Euler_Pressure2D(Q,gas_gamma) <= 0)
    assert(isempty(find(Euler_Pressure2D(Q,gas_gamma) <= 0)));
    
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
    
    figure(5)
    perc_tcells2 = length(ind_neg)*100/K;
    maxp_tcells2 = max(maxp_tcells2,perc_tcells2);
    plot(time,perc_tcells2,'r.')
    xlim([0,FinalTime])
    ylim([0,100])
    xlabel('time')
    ylabel('Percentage of cells flagged')
    title(['Max = ',num2str(maxp_tcells2),'%',' of ',num2str(K),' cells'])
    hold all

    
    % Render every 25 time steps
    if(~mod(tstep,5) || time == FinalTime)
        
        fprintf('   ---> Time = %f, CFL = %f, dt = %f, number of time steps = %f\n',time,CFL,dt,tstep);
        
        figure(2);
        plot(xavg(ind), yavg(ind), 'r.');
        hold all
        plot(xavg(ind_neg), yavg(ind_neg), 'bx');
        xlim([min(min(x)),max(max(x))])
        ylim([min(min(y)),max(max(y))])
        title(['Troubled cells t=',num2str(time)]);
        drawnow;
        hold off
        
        figure(3)
        PlotField2D(Q(:,:,1),interp,TRI,xout,yout); axis tight; drawnow;
        
        figure(4)
        PlotField2D(Euler_Pressure2D(Q,gas_gamma),interp,TRI,xout,yout); axis tight; drawnow;
        
        
        
        
        pause(0.1);
        
    end
    
    if(Save_soln)
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
    
    if(isempty(fixed_dt))
        dt = EulerDT2D(Q,gas_gamma);
    else
        dt = fixed_dt;
        CFL = EulerCFL2D(Q,dt,gas_gamma);
    end
    
    tstep = tstep+1;
end

if(Save_soln)
    Q_save{1,stime}   = Q;
    ind_save{1,stime} = ind;
end

return;
