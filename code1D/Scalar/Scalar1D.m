function [u] = Scalar1D(u,ind_fname)

% Purpose  : Integrate 1D Scalar equation until 
%            FinalTime starting with
%            initial condition u in the domain [xL,xR].

Globals1D_DG;
Globals1D_MLP;

% Set exact flux and Jacobian
[flux,dflux] = Set_scalar_flux1D(model);

time = 0;

xmin = min(abs(x(1,:)-x(2,:)));
iter = 0;
xcen = mean(x,1);

if(save_ind)
    fid = fopen(ind_fname,'w');
end



% Limit initial solution
ind0  = Scalar1D_Tcells(u,bc_cond);
u     = Scalar1D_limit(u,ind0,bc_cond);
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



% outer time step loop 
while(time<FinalTime)
  
  speed = max(max(abs(dflux(u))));
  dt = CFL* min(xmin/abs(speed));
  
  if(time+dt>FinalTime)
    dt = FinalTime-time;
  end

  % 3rd order SSP Runge-Kutta
  
  % SSP RK Stage 1.
  rhsu  = ScalarRHS1D_weak(u,flux,dflux);
  u1  = u  + dt*rhsu;

  % Limit fields
  ind1  = Scalar1D_Tcells(u1,bc_cond);
  u1    = Scalar1D_limit(u1,ind1,bc_cond);
  if(save_ind)
      Tcell_write1D(fid,time+dt,xcen(ind1));
  end
  

  % SSP RK Stage 2.
  rhsu  = ScalarRHS1D_weak(u1,flux,dflux);
  u2   = (3*u  + u1  + dt*rhsu )/4;

  % Limit fields
  ind2  = Scalar1D_Tcells(u2,bc_cond);
  u2    = Scalar1D_limit(u2,ind2,bc_cond);
  if(save_ind)
      Tcell_write1D(fid,time+dt,xcen(ind2));
  end

  % SSP RK Stage 3.
  rhsu  = ScalarRHS1D_weak(u2,flux,dflux);
  u  = (u  + 2*u2  + 2*dt*rhsu )/3;
  
  % Limit solution
  ind3  = Scalar1D_Tcells(u,bc_cond);
  u     = Scalar1D_limit(u,ind3,bc_cond);
  if(save_ind)
      Tcell_write1D(fid,time+dt,xcen(ind3));
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
      figure(1)
      plot(x,u,'b-','LineWidth',2)
      hold all
      xlabel('x')
      title(['time = ',num2str(time)])
      hold off

      pause(0.1)
  end
    
  iter = iter + 1;
  
end

if(save_ind)
    fclose(fid);
end

return
