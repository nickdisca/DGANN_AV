function [u] = Scalar1D(u,Problem,Mesh,Limit,Net,Viscosity,NetVisc,Output)

% Purpose  : Integrate 1D Scalar equation until
%            FinalTime starting with
%            initial condition u in the domain [xL,xR].

% Set exact flux and Jacobian
[flux,dflux] = Set_scalar_flux1D(Problem.model);

time = 0;
dt = 0;

xmin = min(abs(Mesh.x(1,:)-Mesh.x(2,:)));
iter = 0;
xcen = mean(Mesh.x,1);

if(Output.save_ind)
    fid = fopen(strcat(Output.fname_base,'_tcells.dat'),'w');
end
if(Output.save_visc)
    fid2 = fopen(strcat(Output.fname_base,'_visc.dat'),'w');
end


% Limit initial solution
ind0  = Scalar1D_Tcells(u,Problem.bc_cond,Mesh,Limit,Net);
u     = Scalar1D_limit(u,ind0,Problem.bc_cond,Limit.Limiter,Mesh);
if(Output.save_ind)
    Tcell_write1D(fid,time,xcen(ind0));
end

% Initialize solution at previous time step
u_tmp=zeros(length(u(:)),3);


% outer time step loop
while(time<Problem.FinalTime)
    
    %Compute artificial viscosity
    uold = reshape(u_tmp(:,1),Mesh.Np,Mesh.K);
    mu_piece = Scalar1D_viscosity(u, uold, dflux, Viscosity, Problem, Mesh, dt, iter, NetVisc);
    mu_piece=max(mu_piece,0);
    mu_vals=Scalar1D_smooth_viscosity(mu_piece,Mesh.x); 
    mu_vals=max(mu_vals,0);
    maxvisc=max(abs(mu_vals(:)));
    
    if(Output.save_visc && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
        Visc_write1D(fid2,time,mu_vals);
    end
    
    %Set timestep
    speed = max(max(abs(dflux(u))));
    dt = Problem.CFL*1/(speed*Mesh.N^2/min(Mesh.hK)+maxvisc*Mesh.N^4/min(Mesh.hK)^2);
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % Save aux variable (needed for EV)
    u_tmp(:,2)=u(:);
    
    % 3rd order SSP Runge-Kutta
    if strcmp(Problem.RK,'SSP3')
        
        rhsu  = ScalarRHS1D_weak(u,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
        u1  = u  + dt*rhsu;
        
        % Limit fields
        ind1  = Scalar1D_Tcells(u1,Problem.bc_cond,Mesh,Limit,Net);
        u1    = Scalar1D_limit(u1,ind1,Problem.bc_cond,Limit.Limiter,Mesh);
        if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
            Tcell_write1D(fid,time+dt,xcen(ind1));
        end 
        
        % SSP RK Stage 2.
        rhsu  = ScalarRHS1D_weak(u1,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
        u2   = (3*u  + u1  + dt*rhsu )/4;
        
        % Limit fields
        ind2  = Scalar1D_Tcells(u2,Problem.bc_cond,Mesh,Limit,Net);
        u2    = Scalar1D_limit(u2,ind2,Problem.bc_cond,Limit.Limiter,Mesh);
        if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
            Tcell_write1D(fid,time+dt,xcen(ind2));
        end
        
        % SSP RK Stage 3.
        rhsu  = ScalarRHS1D_weak(u2,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
        u  = (u  + 2*u2  + 2*dt*rhsu )/3;
        
        % Limit solution
        ind3  = Scalar1D_Tcells(u,Problem.bc_cond,Mesh,Limit,Net);
        u     = Scalar1D_limit(u,ind3,Problem.bc_cond,Limit.Limiter,Mesh);
        if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
            Tcell_write1D(fid,time+dt,xcen(ind3));
        end
       
        
        
    % 4th order low storage Runge-Kutta    
    elseif strcmp(Problem.RK,'LS54')
        
        A_RK(1)=0; A_RK(2)=-0.4178904745; A_RK(3)=-1.192151694643; A_RK(4)=-1.697784692471; A_RK(5)=-1.514183444257;
        B_RK(1)=0.1496590219993; B_RK(2)=0.3792103129999; B_RK(3)=0.8229550293869; B_RK(4)=0.6994504559488; B_RK(5)=0.1530572479681;
        c_RK(1)=0; c_RK(2)=0.1496590219993; c_RK(3)=0.3704009573644; c_RK(4)=0.6222557631345; c_RK(5)=0.9582821306748;
        nsteps=5;
        
        UU=u; VV=zeros(size(u));
        for index=1:nsteps
            
            rhsu  = ScalarRHS1D_weak(UU,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
            
            VV=A_RK(index)*VV+dt*rhsu;
            
            UU=UU+B_RK(index)*VV;
            
            ind_substage  = Scalar1D_Tcells(UU,Problem.bc_cond,Mesh,Limit,Net);
            
            UU  = Scalar1D_limit(UU,ind_substage,Problem.bc_cond,Limit.Limiter,Mesh);
            if(Output.save_ind && (mod(iter,Output.save_iter) == 0 || time+dt >= Problem.FinalTime))
                Tcell_write1D(fid,time+dt,xcen(ind_substage));
            end
        end
            
        u=UU;
        
        
        
    else
        
        error('Time integration scheme not defined');
    end
    
    % Increment time and iteration
    time = time+dt;
    iter = iter+1;
    
    % Increment saved variables (needed for EV)
    u_tmp(:,3)=u(:); u_tmp(:,1)=u_tmp(:,2); u_tmp(:,2)=u_tmp(:,3);
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        figure(1)
        plot(Mesh.x(:),u(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('u')
        title(['time = ',num2str(time)])
        
        pause(0.1)
    end
    
        
end

if(Output.save_ind)
    fclose(fid);
end

return
