function q = Euler1D(q,gas_gamma, gas_const,ind_fname)

% Purpose  : Integrate 1D Euler equations until FinalTime

Globals1D_DG;
Globals1D_MLP;

time = 0;

xmin = min(abs(x(1,:)-x(2,:)));
iter = 0;
xcen = mean(x,1);

% Limit initial solution
if(save_ind)
    fid = fopen(ind_fname,'w');
end
ind0 = Tcells_Euler_type(q,gas_gamma, gas_const,ind_var);
q    = SlopeLimit_Euler_type(q,gas_gamma, gas_const,lim_var,ind0);

if(save_ind)
    Tcell_write(fid,time,xcen(ind0));
    figure(100)
    subplot(1,3,1)
    plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
    xlabel('x')
    ylabel('t')
    xlim([bnd_l bnd_r])
    ylim([0 FinalTime])
    title('RK stage 1')
    hold all
    subplot(1,3,2)
    plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
    xlabel('x')
    ylabel('t')
    xlim([bnd_l bnd_r])
    ylim([0 FinalTime])
    title('RK stage 2')
    hold all
    subplot(1,3,3)
    plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
    xlabel('x')
    ylabel('t')
    xlim([bnd_l bnd_r])
    ylim([0 FinalTime])
    title('RK stage 2')
    hold all
end

iter = 0;

% outer time step loop
while(time<FinalTime)
    
    pre    = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
    lambda = sqrt(gas_gamma*pre./q(:,:,1)) + abs(q(:,:,2)./q(:,:,1));
    dt     = CFL*min(min(xmin./(lambda)));
    
    if(time+dt>FinalTime)
        dt = FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    [rhsq]  = EulerRHS1D_weak(q, gas_gamma, gas_const);
    q1      = q + dt*rhsq;
    
    ind1 = Tcells_Euler_type(q1,gas_gamma, gas_const,ind_var);
    q1   = SlopeLimit_Euler_type(q1,gas_gamma, gas_const,lim_var,ind1);
   
    
    if(save_ind)
        Tcell_write(fid,time+dt,xcen(ind1));
    end
    
    pre = (gas_gamma-1)*(q1(:,:,3) - 0.5*q1(:,:,2).^2./q1(:,:,1));
    if( min(min(real(q1(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 2.
    [rhsq]  = EulerRHS1D_weak(q1, gas_gamma, gas_const);
    q2      = (3*q + (q1 + dt*rhsq))/4.0;
    
    ind2 = Tcells_Euler_type(q2,gas_gamma, gas_const,ind_var);
    q2   = SlopeLimit_Euler_type(q2,gas_gamma, gas_const,lim_var,ind2);
    
    
    if(save_ind)
        Tcell_write(fid,time+dt,xcen(ind2));
    end
    
    pre    = (gas_gamma-1)*(q2(:,:,3) - 0.5*q2(:,:,2).^2./q2(:,:,1));
    if( min(min(real(q2(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 3.
    [rhsq]  = EulerRHS1D_weak(q2, gas_gamma, gas_const);
    q       = (q + 2*(q2 + dt*rhsq))/3.0;
    
    ind3 = Tcells_Euler_type(q,gas_gamma, gas_const,ind_var);
    q    = SlopeLimit_Euler_type(q,gas_gamma, gas_const,lim_var,ind3);    

    if(save_ind)
        Tcell_write(fid,time+dt,xcen(ind3));
    end

 
    pre    = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
    if( min(min(real(q(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    % Increment time and adapt timestep
    time = time+dt;
    
    if(save_ind)
        figure(100)
        subplot(1,3,1)
        plot(xcen(ind1),ones(1,length(ind1))*time,'r.')
        subplot(1,3,2)
        plot(xcen(ind2),ones(1,length(ind2))*time,'r.')
        subplot(1,3,3)
        plot(xcen(ind3),ones(1,length(ind3))*time,'r.')
    end
    
    
    if(mod(iter,plot_iter) == 0 || time >= FinalTime)
        depth = q(:,:,1);
        vel   = q(:,:,2)./q(:,:,1);
        figure(1)
        subplot(3,1,1)
        plot(x(:),depth(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Density')
        title(['time = ',num2str(time)])
        
        subplot(3,1,2)
        plot(x(:),vel(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Velocity')

        subplot(3,1,3)
        plot(x(:),pre(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Pressure')

        pause(.1)
        
    end
    
    
    iter = iter + 1;
    
end

if(save_ind)
    fclose(fid);
end

return
