function PlotEuler1D(fname_base,x_ran,var_ran,t_ran,ref_avail,ref_fname,x,RKscheme)


%plot solution at final time
soln = load(strcat(fname_base,'.dat'));

if(ref_avail) 
    ref_soln = load(ref_fname);
end

figure(1); clf;
figure(2); clf;
figure(3); clf;

if(ref_avail)
    figure(1)
    plot(ref_soln(:,1),ref_soln(:,2),'k-','LineWidth',2) 
    hold all
    
    figure(2)
    plot(ref_soln(:,1),ref_soln(:,3),'k-','LineWidth',2) 
    hold all
    
    figure(3)
    plot(ref_soln(:,1),ref_soln(:,4),'k-','LineWidth',2) 
    hold all
end

figure(1)
plot(soln(:,1),soln(:,2),'r-','LineWidth',2)
figure(2)
plot(soln(:,1),soln(:,3),'r-','LineWidth',2)
figure(3)
plot(soln(:,1),soln(:,4),'r-','LineWidth',2)

if(ref_avail)
    figure(1)
    legend('Reference','Numerical')
    
    figure(2)
    legend('Reference','Numerical')
    
    figure(3)
    legend('Reference','Numerical')
end

figure(1)
xlim(x_ran);
ylim(var_ran(1,:));
xlabel({'$x$'},'interpreter','latex')
ylabel({'$\rho$'},'interpreter','latex')
set(gca,'FontSize',20)
fname = sprintf('%s_soln_density.pdf',fname_base);
print(fname,'-dpdf')

figure(2)
xlim(x_ran);
ylim(var_ran(2,:));
xlabel({'$x$'},'interpreter','latex')
ylabel({'$v$'},'interpreter','latex')
set(gca,'FontSize',20)
fname = sprintf('%s_soln_vel.pdf',fname_base);
print(fname,'-dpdf')

figure(3)
xlim(x_ran);
ylim(var_ran(3,:));
xlabel({'$x$'},'interpreter','latex')
ylabel({'$p$'},'interpreter','latex')
set(gca,'FontSize',20)
fname = sprintf('%s_soln_pre.pdf',fname_base);
print(fname,'-dpdf')


fid = fopen(strcat(fname_base,'_tcells.dat'));

tline = fgetl(fid);
data  = str2num(char(regexp(tline,', ','split')));
t     = data(1);
ind   = data(2:end);

figure;
clf
plot(t*ones(1,length(ind)),ind,'r.')
xlim(x_ran)
ylim(t_ran)
xlabel({'$x$'},'interpreter','latex')
ylabel({'$t$'},'interpreter','latex')
set(gca,'FontSize',20)
hold all

tline = fgetl(fid);
while ischar(tline)
    
    ind_i=[];
    
    for i=1:(3*strcmp(RKscheme,'SSP3')+5*strcmp(RKscheme,'LS54'))
        data  = str2num(char(regexp(tline,', ','split')));
        t     = data(1);
        ind_i = [ind_i; data(2:end)];
        tline = fgetl(fid);
    end
    
    inda = unique(ind_i(:));
    plot(inda,t*ones(1,length(inda)),'r.');
    
end

fname = sprintf('%s_tcells.pdf',fname_base);
print(fname,'-dpdf')


%plot viscosity
figure;
fid = fopen(strcat(fname_base,'_visc.dat'));

tline = fgetl(fid);
visc_= [];
t = [];
while ischar(tline)

    data  = str2num(char(regexp(tline,', ','split')));
    t     = [t; data(1)];
    visc_ = [visc_, data(2:end)];
    tline = fgetl(fid);
    
end

visc_=log(max(visc_,1e-5))/log(10);
plot_visc=pcolor(x(:),t,visc_'); shading interp; set(plot_visc, 'EdgeColor', 'none'); colorbar; colormap jet;
xlabel({'$x$'},'interpreter','latex','FontSize',20); 
ylabel({'$t$'},'interpreter','latex','FontSize',20);

fname = sprintf('%s_visc.pdf',fname_base);
print(fname,'-dpdf')




fclose(fid);



