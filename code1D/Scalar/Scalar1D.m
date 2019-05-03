function [u] = Scalar1D(u,Problem,Mesh,Limit,Net,Output)

% Purpose  : Integrate 1D Scalar equation until
%            FinalTime starting with
%            initial condition u in the domain [xL,xR].

% Globals1D_DG;
% Globals1D_MLP;

% Set exact flux and Jacobian
[flux,dflux] = Set_scalar_flux1D(Problem.model);

time = 0;

xmin = min(abs(Mesh.x(1,:)-Mesh.x(2,:)));
iter = 0;
xcen = mean(Mesh.x,1);

if(Output.save_ind)
    fid = fopen(strcat(Output.fname_base,'_tcells.dat'),'w');
end



% Limit initial solution
ind0  = Scalar1D_Tcells(u,Problem.bc_cond,Mesh,Limit,Net);
u     = Scalar1D_limit(u,ind0,Problem.bc_cond,Limit.Limiter,Mesh);
if(Output.save_ind)
    Tcell_write1D(fid,time,xcen(ind0));
    figure(10)
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


% iter_p = 0;
% figure(100)
% plot(x(:),u(:),'b-','LineWidth',2)
% xlabel('x')
% ylabel('u')
% title(['time = ',num2str(time)])
% ylim([-1.5,3.5])
% set(gca,'FontSize',20)
% print(['soln_',num2str(iter_p),'.pdf'],'-dpdf')
% 
% figure(200)
% ind_all = ind0;
% plot(xcen(ind_all),ones(1,length(ind_all))*time,'b.')
% xlabel('x')
% ylabel('t')
% xlim([bnd_l bnd_r])
% ylim([0 FinalTime])
% title('Cells Flagged')
% set(gca,'FontSize',20)
% hold all
% print(['tcells_',num2str(iter_p),'.pdf'],'-dpdf')

% figure(200)
% ind_all = ind0;
% plot(xcen(ind_all),ones(1,length(ind_all))*time,'b.')
% xlabel('x')
% ylabel('t')
% xlim([bnd_l bnd_r])
% ylim([0 FinalTime])
% title('Cells Flagged')
% set(gca,'FontSize',20)
% hold all
% print(['tcells_',num2str(iter_p),'.pdf'],'-dpdf')


% outer time step loop
while(time<Problem.FinalTime)
    
    speed = max(max(abs(dflux(u))));
    dt = Problem.CFL* min(xmin/abs(speed));
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    rhsu  = ScalarRHS1D_weak(u,flux,dflux,Problem.bc_cond,Mesh);
    u1  = u  + dt*rhsu;
    
    % Limit fields
    ind1  = Scalar1D_Tcells(u1,Problem.bc_cond,Mesh,Limit,Net);
    u1    = Scalar1D_limit(u1,ind1,Problem.bc_cond,Limit.Limiter,Mesh);
    if(Output.save_ind && (mod(iter,Output.plot_iter) == 0 || time+dt >= Problem.FinalTime))
        Tcell_write1D(fid,time+dt,xcen(ind1));
    end
    
    
    % SSP RK Stage 2.
    rhsu  = ScalarRHS1D_weak(u1,flux,dflux,Problem.bc_cond,Mesh);
    u2   = (3*u  + u1  + dt*rhsu )/4;
    
    % Limit fields
    ind2  = Scalar1D_Tcells(u2,Problem.bc_cond,Mesh,Limit,Net);
    u2    = Scalar1D_limit(u2,ind2,Problem.bc_cond,Limit.Limiter,Mesh);
    if(Output.save_ind && (mod(iter,Output.plot_iter) == 0 || time+dt >= Problem.FinalTime))
        Tcell_write1D(fid,time+dt,xcen(ind2));
    end
    
    % SSP RK Stage 3.
    rhsu  = ScalarRHS1D_weak(u2,flux,dflux,Problem.bc_cond,Mesh);
    u  = (u  + 2*u2  + 2*dt*rhsu )/3;
    
    % Limit solution
    ind3  = Scalar1D_Tcells(u,Problem.bc_cond,Mesh,Limit,Net);
    u     = Scalar1D_limit(u,ind3,Problem.bc_cond,Limit.Limiter,Mesh);
    if(Output.save_ind && (mod(iter,Output.plot_iter) == 0 || time+dt >= Problem.FinalTime))
        Tcell_write1D(fid,time+dt,xcen(ind3));
    end
    
    % Increment time and adapt timestep
    time = time+dt;
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        figure(1)
        plot(Mesh.x(:),u(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('u')
        title(['time = ',num2str(time)])
        
        figure(10)
        subplot(1,3,1)
        plot(xcen(ind1),ones(1,length(ind1))*time,'r.')
        subplot(1,3,2)
        plot(xcen(ind2),ones(1,length(ind2))*time,'r.')
        subplot(1,3,3)
        plot(xcen(ind3),ones(1,length(ind3))*time,'r.')
        
        pause(0.1)
    end
    
%     if(mod(iter,5) == 0 || time >= FinalTime)
%         iter_p = iter_p + 1;
%         figure(100)
%         plot(x(:),u(:),'b-','LineWidth',2)
%         xlabel('x')
%         ylabel('u')
%         title(['time = ',num2str(time)])
%         ylim([-1.5,3.5])
%         set(gca,'FontSize',20)
%         print(['soln_',num2str(iter_p),'.pdf'],'-dpdf')
%         
%         figure(200)
%         ind_all = unique([ind1,ind2,ind3]);
%         plot(xcen(ind_all),ones(1,length(ind_all))*time,'b.')
%         xlabel('x')
%         ylabel('t')
%         xlim([bnd_l bnd_r])
%         ylim([0 FinalTime])
%         title('Cells Flagged')
%         set(gca,'FontSize',20)
%         hold all
%         print(['tcells_',num2str(iter_p),'.pdf'],'-dpdf')
%         
%         pause(.1)
%     end
    
    iter = iter + 1;
    
end

if(Output.save_ind)
    fclose(fid);
end

return
