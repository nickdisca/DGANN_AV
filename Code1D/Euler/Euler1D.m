function q = Euler1D(q,Problem,Mesh,Limit,Net,Viscosity,NetVisc,Output)

% Purpose  : Integrate 1D Euler equations until FinalTime

time = 0;
dt = 0;

xmin = min(abs(Mesh.x(1,:)-Mesh.x(2,:)));
iter = 0;
xcen = mean(Mesh.x,1);

% Limit initial solution
if(Output.save_ind)
    fid = fopen(strcat(Output.fname_base,'_tcells.dat'),'w');
end
if(Output.save_visc)
    fid2 = fopen(strcat(Output.fname_base,'_visc.dat'),'w');
end

ind0 = Tcells_Euler_type1D(q,Problem,Mesh,Limit,Net);
q    = SlopeLimit_Euler_type1D(q,ind0,Problem,Limit,Mesh);

if(Output.save_ind)
    Tcell_write1D(fid,time,xcen(ind0));
end

% Initialize solution at previous time step
q_tmp=zeros(length(q(:)),3);


% outer time step loop
while(time<Problem.FinalTime)
    
    %Compute artificial viscosity
    qold = reshape(q_tmp(:,1),Mesh.Np,Mesh.K,3);
    mu_piece = Euler1D_viscosity(q, qold, Viscosity, Problem, Mesh, dt, iter, NetVisc);
    mu_piece=max(mu_piece,0);
    mu_vals=Scalar1D_smooth_viscosity(mu_piece,Mesh.x); 
    mu_vals=max(mu_vals,0);
    maxvisc=max(abs(mu_vals(:)));
    
    if(Output.save_visc && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
        Visc_write1D(fid2,time,mu_vals);
    end
    
    %Apply same AV for all equations
    mu_vals=repmat(mu_vals,1,1,3);
    
    %Set timestep
    pre    = (Problem.gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
    c_sound = sqrt(Problem.gas_gamma*pre./q(:,:,1));
    %c_sound = abs(c_sound);
    if min(min(q(:,:,1)))<=0 || min(pre(:))<=0
        warning('Positivity loss');
    end
    lambda = c_sound + abs(q(:,:,2)./q(:,:,1));
    dt = Problem.CFL*1/(max(lambda(:))*Mesh.N^2/min(Mesh.hK)+maxvisc*Mesh.N^4/min(Mesh.hK)^2);
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % Save aux variable (needed for EV)
    q_tmp(:,2)=q(:);
    
    % 3rd order SSP Runge-Kutta
    if strcmp(Problem.RK,'SSP3')
    
        % SSP RK Stage 1.
        [rhsq]  = EulerRHS1D_weak(q, Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
        q1      = q + dt*rhsq;
        
        %Limit fields
        ind1 = Tcells_Euler_type1D(q1,Problem,Mesh,Limit,Net);
        q1   = SlopeLimit_Euler_type1D(q1,ind1,Problem,Limit,Mesh);
        if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
            Tcell_write1D(fid,time+dt,xcen(ind1));
        end
        
        
        % SSP RK Stage 2.
        [rhsq]  = EulerRHS1D_weak(q1, Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
        q2      = (3*q + (q1 + dt*rhsq))/4.0;
        
        %Limit fields
        ind2 = Tcells_Euler_type1D(q2,Problem,Mesh,Limit,Net);
        q2   = SlopeLimit_Euler_type1D(q2,ind2,Problem,Limit,Mesh);
        if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
            Tcell_write1D(fid,time+dt,xcen(ind2));
        end
        
        
        % SSP RK Stage 3.
        [rhsq]  = EulerRHS1D_weak(q2,Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
        q       = (q + 2*(q2 + dt*rhsq))/3.0;
        
        %Limit fields
        ind3 = Tcells_Euler_type1D(q,Problem,Mesh,Limit,Net);
        q    = SlopeLimit_Euler_type1D(q,ind3,Problem,Limit,Mesh);
        if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
            Tcell_write1D(fid,time+dt,xcen(ind3));
        end
        
    
    % 4th order low storage Runge-Kutta    
    elseif strcmp(Problem.RK,'LS54')
        
        A_RK(1)=0; A_RK(2)=-0.4178904745; A_RK(3)=-1.192151694643; A_RK(4)=-1.697784692471; A_RK(5)=-1.514183444257;
        B_RK(1)=0.1496590219993; B_RK(2)=0.3792103129999; B_RK(3)=0.8229550293869; B_RK(4)=0.6994504559488; B_RK(5)=0.1530572479681;
        c_RK(1)=0; c_RK(2)=0.1496590219993; c_RK(3)=0.3704009573644; c_RK(4)=0.6222557631345; c_RK(5)=0.9582821306748;
        nsteps=5;
        
        QQ=q; VV=zeros(size(q));
        for index=1:nsteps
            
            rhsq  = EulerRHS1D_weak(QQ,Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
            
            VV=A_RK(index)*VV+dt*rhsq;
            
            QQ=QQ+B_RK(index)*VV;
            
            ind_substage  = Tcells_Euler_type1D(QQ,Problem,Mesh,Limit,Net);
            
            QQ  =  SlopeLimit_Euler_type1D(QQ,ind_substage,Problem,Limit,Mesh);
            if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
                Tcell_write1D(fid,time+dt,xcen(ind_substage));
            end
        end
            
        q=QQ;
        
    else
        
        error('Time integration scheme not defined');
    end
    
    % Increment time and adapt timestep
    time = time+dt;   
    iter = iter + 1;
    
    % Increment saved variables (needed for EV)
    q_tmp(:,3)=q(:); q_tmp(:,1)=q_tmp(:,2); q_tmp(:,2)=q_tmp(:,3);
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        density = q(:,:,1);
        vel   = q(:,:,2)./q(:,:,1);
        figure(1)
        subplot(3,1,1)
        plot(Mesh.x(:),density(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Density')
        title(['time = ',num2str(time)])
        
        subplot(3,1,2)
        plot(Mesh.x(:),vel(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Velocity')
        
        subplot(3,1,3)
        plot(Mesh.x(:),pre(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Pressure')
        
        pause(.1)
        
    end
    

      
end

if(Output.save_ind)
    fclose(fid);
end

return
