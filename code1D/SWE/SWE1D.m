function q = SWE1D(q,gravity,ind_fname)

% Purpose  : Integrate 1D Shallow water equations until FinalTime

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
ind0 = Tcells_SWE_type1D(q,ind_var,bc_cond);
q    = SlopeLimit_SWE_type1D(q,gravity,lim_var,ind0,bc_cond);

% [qc,Char] = SWEUtoC1D(q,gravity);
% [Char(:,:,1),ind0x] = SlopeLimitN_SWE2(Char(:,:,1));
% [Char(:,:,2),ind0y] = SlopeLimitN_SWE2(Char(:,:,2));
% [q] = SWECtoU1D(Char,qc,gravity);
% ind0 = unique([ind0x,ind0y]);
if(save_ind)
    Tcell_write1D(fid,time,xcen(ind0));
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
    
    lambda = sqrt(gravity*q(:,:,1)) + abs(q(:,:,2)./q(:,:,1));
    dt = CFL*min(min(xmin./(lambda)));
    
    if(time+dt>FinalTime)
        dt = FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    [rhsq]  = SWERHS1D_weak(q, gravity);
    q1      = q + dt*rhsq;
    
    ind1 = Tcells_SWE_type1D(q1,ind_var,bc_cond);
    q1   = SlopeLimit_SWE_type1D(q1,gravity,lim_var,ind1,bc_cond);
    
%     [qc,Char] = SWEUtoC1D(q1,gravity);
%     [Char(:,:,1),ind1x] = SlopeLimitN_SWE2(Char(:,:,1));
%     [Char(:,:,2),ind1y] = SlopeLimitN_SWE2(Char(:,:,2));
%     [q1] = SWECtoU1D(Char,qc,gravity);
%     ind1 = unique([ind1x,ind1y]);
    
    if(save_ind)
        Tcell_write1D(fid,time+dt,xcen(ind1));
    end
    
    if( min(min(real(q1(:,:,1)))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 2.
    [rhsq]  = SWERHS1D_weak(q1, gravity);
    q2      = (3*q + (q1 + dt*rhsq))/4.0;
    
    ind2 = Tcells_SWE_type1D(q2,ind_var,bc_cond);
    q2   = SlopeLimit_SWE_type1D(q2,gravity,lim_var,ind2,bc_cond);
    
%     [qc,Char] = SWEUtoC1D(q2,gravity);
%     [Char(:,:,1),ind2x] = SlopeLimitN_SWE2(Char(:,:,1));
%     [Char(:,:,2),ind2y] = SlopeLimitN_SWE2(Char(:,:,2));
%     [q2] = SWECtoU1D(Char,qc,gravity);
%     ind2 = unique([ind2x,ind2y]);
    
    if(save_ind)
        Tcell_write1D(fid,time+dt,xcen(ind2));
    end
    
    if( min(min(real(q2(:,:,1)))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 3.
    [rhsq]  = SWERHS1D_weak(q2, gravity);
    q       = (q + 2*(q2 + dt*rhsq))/3.0;
    
    ind3 = Tcells_SWE_type1D(q,ind_var,bc_cond);
    q    = SlopeLimit_SWE_type1D(q,gravity,lim_var,ind3,bc_cond);

%     [qc,Char] = SWEUtoC1D(q1,gravity);
%     [Char(:,:,1),ind3x] = SlopeLimitN_SWE2(Char(:,:,1));
%     [Char(:,:,2),ind3y] = SlopeLimitN_SWE2(Char(:,:,2));
%     [q] = SWECtoU1D(Char,qc,gravity);
%     ind3 = unique([ind3x,ind3y]);    

    if(save_ind)
        Tcell_write1D(fid,time+dt,xcen(ind3));
    end
    
    if( min(min(real(q(:,:,1)))) <= 0.0)
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
        subplot(2,1,1)
        plot(x(:),depth(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Depth')
        title(['time = ',num2str(time)])
        
        subplot(2,1,2)
        plot(x(:),vel(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Velocity')

        pause(.1)
        
    end
    
    
    iter = iter + 1;
    
end

if(save_ind)
    fclose(fid);
end

return
