function q = SWE1D(q,Problem,Mesh,Limit,Net,Output)

% Purpose  : Integrate 1D Shallow water equations until FinalTime

% Globals1D_DG;
% Globals1D_MLP;

time = 0;

xmin = min(abs(Mesh.x(1,:)-Mesh.x(2,:)));
iter = 0;
xcen = mean(Mesh.x,1);

% Limit initial solution
if(Output.save_ind)
    fid = fopen(strcat(Output.fname_base,'_tcells.dat'),'w');
end

ind0 = Tcells_SWE_type1D(q,Problem.bc_cond,Mesh,Limit,Net);
q    = SlopeLimit_SWE_type1D(q,ind0,Problem.gravity,Problem.bc_cond,Limit,Mesh);

% [qc,Char] = SWEUtoC1D(q,gravity);
% [Char(:,:,1),ind0x] = SlopeLimitN_SWE2(Char(:,:,1));
% [Char(:,:,2),ind0y] = SlopeLimitN_SWE2(Char(:,:,2));
% [q] = SWECtoU1D(Char,qc,gravity);
% ind0 = unique([ind0x,ind0y]);
if(Output.save_ind)
    Tcell_write1D(fid,time,xcen(ind0));
    figure(100)
    subplot(1,3,1)
    plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
    xlabel('x')
    ylabel('t')
    xlim([Mesh.bnd_l Mesh.bnd_r])
    ylim([0 Problem.FinalTime])
    title('RK stage 1')
    hold all
    subplot(1,3,2)
    plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
    xlabel('x')
    ylabel('t')
    xlim([Mesh.bnd_l Mesh.bnd_r])
    ylim([0 Problem.FinalTime])
    title('RK stage 2')
    hold all
    subplot(1,3,3)
    plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
    xlabel('x')
    ylabel('t')
    xlim([Mesh.bnd_l Mesh.bnd_r])
    ylim([0 Problem.FinalTime])
    title('RK stage 3')
    hold all
end


% outer time step loop
while(time<Problem.FinalTime)
    
    lambda = sqrt(Problem.gravity*q(:,:,1)) + abs(q(:,:,2)./q(:,:,1));
    dt = Problem.CFL*min(min(xmin./(lambda)));
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    [rhsq]  = SWERHS1D_weak(q,Problem.gravity,Problem.bc_cond,Mesh);
    q1      = q + dt*rhsq;
    
    ind1 = Tcells_SWE_type1D(q1,Problem.bc_cond,Mesh,Limit,Net);
    q1   = SlopeLimit_SWE_type1D(q1,ind1,Problem.gravity,Problem.bc_cond,Limit,Mesh);
    
%     [qc,Char] = SWEUtoC1D(q1,gravity);
%     [Char(:,:,1),ind1x] = SlopeLimitN_SWE2(Char(:,:,1));
%     [Char(:,:,2),ind1y] = SlopeLimitN_SWE2(Char(:,:,2));
%     [q1] = SWECtoU1D(Char,qc,gravity);
%     ind1 = unique([ind1x,ind1y]);
    
    if(Output.save_ind && (mod(iter,Output.plot_iter) == 0 || time+dt >= Problem.FinalTime))
        Tcell_write1D(fid,time+dt,xcen(ind1));
    end
    
    if( min(min(real(q1(:,:,1)))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 2.
    [rhsq]  = SWERHS1D_weak(q1,Problem.gravity,Problem.bc_cond,Mesh);
    q2      = (3*q + (q1 + dt*rhsq))/4.0;
    
    ind2 = Tcells_SWE_type1D(q2,Problem.bc_cond,Mesh,Limit,Net);
    q2   = SlopeLimit_SWE_type1D(q2,ind2,Problem.gravity,Problem.bc_cond,Limit,Mesh);
    
%     [qc,Char] = SWEUtoC1D(q2,gravity);
%     [Char(:,:,1),ind2x] = SlopeLimitN_SWE2(Char(:,:,1));
%     [Char(:,:,2),ind2y] = SlopeLimitN_SWE2(Char(:,:,2));
%     [q2] = SWECtoU1D(Char,qc,gravity);
%     ind2 = unique([ind2x,ind2y]);
    
    if(Output.save_ind && (mod(iter,Output.plot_iter) == 0 || time+dt >= Problem.FinalTime))
        Tcell_write1D(fid,time+dt,xcen(ind2));
    end
    
    if( min(min(real(q2(:,:,1)))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 3.
    [rhsq]  = SWERHS1D_weak(q2,Problem.gravity,Problem.bc_cond,Mesh);
    q       = (q + 2*(q2 + dt*rhsq))/3.0;
    
    ind3 = Tcells_SWE_type1D(q,Problem.bc_cond,Mesh,Limit,Net);
    q    = SlopeLimit_SWE_type1D(q,ind3,Problem.gravity,Problem.bc_cond,Limit,Mesh);

%     [qc,Char] = SWEUtoC1D(q1,gravity);
%     [Char(:,:,1),ind3x] = SlopeLimitN_SWE2(Char(:,:,1));
%     [Char(:,:,2),ind3y] = SlopeLimitN_SWE2(Char(:,:,2));
%     [q] = SWECtoU1D(Char,qc,gravity);
%     ind3 = unique([ind3x,ind3y]);    

    if(Output.save_ind && (mod(iter,Output.plot_iter) == 0 || time+dt >= Problem.FinalTime))
        Tcell_write1D(fid,time+dt,xcen(ind3));
    end
    
    if( min(min(real(q(:,:,1)))) <= 0.0)
        error('Positivity loss!!');
    end
    
    % Increment time and adapt timestep
    time = time+dt;
    
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        depth = q(:,:,1);
        vel   = q(:,:,2)./q(:,:,1);
        figure(1)
        subplot(2,1,1)
        plot(Mesh.x(:),depth(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Depth')
        title(['time = ',num2str(time)])
        
        subplot(2,1,2)
        plot(Mesh.x(:),vel(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Velocity')
        title(['time = ',num2str(time)])
        
        figure(100)
        subplot(1,3,1)
        plot(xcen(ind1),ones(1,length(ind1))*time,'r.')
        subplot(1,3,2)
        plot(xcen(ind2),ones(1,length(ind2))*time,'r.')
        subplot(1,3,3)
        plot(xcen(ind3),ones(1,length(ind3))*time,'r.')

        pause(.1)
        
    end
    
    
    iter = iter + 1;
    
end

if(Output.save_ind)
    fclose(fid);
end

return
