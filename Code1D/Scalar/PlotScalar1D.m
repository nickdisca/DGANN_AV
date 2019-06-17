function PlotScalar1D(fname_base,x_ran,u_ran,t_ran,ref_avail,ref_fname,x,RKscheme)

%plot solution at final time
soln = load(strcat(fname_base,'.dat'));

if(ref_avail) 
    ref_soln = load(ref_fname);
end

figure;
clf
if(ref_avail)
    plot(ref_soln(:,1),ref_soln(:,2),'k-','LineWidth',2) 
    hold all
end
plot(soln(:,1),soln(:,2),'r-','LineWidth',2)
if(ref_avail)
    legend('Reference','Numerical')
end
xlim(x_ran);
ylim(u_ran);
xlabel({'$x$'},'interpreter','latex')
ylabel({'$u$'},'interpreter','latex')
set(gca,'FontSize',20)
fname = sprintf('%s_soln.pdf',fname_base);
print(fname,'-dpdf')


%plot troubled cells
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



