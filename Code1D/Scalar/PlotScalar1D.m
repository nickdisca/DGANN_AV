function PlotScalar1D(fname_base,x_ran,u_ran,t_ran,rk_comb,ref_avail,ref_fname)


soln = load(strcat(fname_base,'.dat'));

if(ref_avail) 
    ref_soln = load(ref_fname);
end

figure(1)
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
xlabel('x')
ylabel('u')
set(gca,'FontSize',20)
fname = sprintf('%s_soln.pdf',fname_base);
print(fname,'-dpdf')


fid = fopen(strcat(fname_base,'_tcells.dat'));

%t=0
tline = fgetl(fid);
data  = str2num(char(regexp(tline,', ','split')));
t     = data(1);
ind   = data(2:end);

figure(2)
clf
if(~rk_comb)
    subplot(1,3,1)
    plot(ind,t*ones(1,length(ind)),'r.')
    xlim(x_ran)
    ylim(t_ran)
    xlabel('x')
    ylabel('t')
    title('RK stage 1')
    set(gca,'FontSize',20)
    hold all
    subplot(1,3,2)
    plot(ind,t*ones(1,length(ind)),'r.')
    xlim(x_ran)
    ylim(t_ran)
    xlabel('x')
    ylabel('t')
    title('RK stage 2')
    set(gca,'FontSize',20)
    hold all
    subplot(1,3,3)
    plot(ind,t*ones(1,length(ind)),'r.')
    xlim(x_ran)
    ylim(t_ran)
    xlabel('x')
    ylabel('t')
    title('RK stage 3')
    set(gca,'FontSize',20)
    hold all
else
    plot(t*ones(1,length(ind)),ind,'r.')
    xlim(x_ran)
    ylim(t_ran)
    xlabel('x')
    ylabel('t')
    set(gca,'FontSize',20)
    hold all
end

tline = fgetl(fid);
while ischar(tline)
    data  = str2num(char(regexp(tline,', ','split')));
    t     = data(1);
    ind1  = data(2:end);
    
    tline = fgetl(fid);
    data  = str2num(char(regexp(tline,', ','split')));
    ind2  = data(2:end);
    
    tline = fgetl(fid);
    data  = str2num(char(regexp(tline,', ','split')));
    ind3  = data(2:end);
 
    
    if(~rk_comb)
        subplot(1,3,1)
        plot(ind1,t*ones(1,length(ind1)),'r.')
        hold all
        subplot(1,3,2)
        plot(ind2,t*ones(1,length(ind2)),'r.')
        hold all
        subplot(1,3,3)
        plot(ind3,t*ones(1,length(ind3)),'r.')
        hold all
    else
        inda = unique([ind1;ind2;ind3]);
        plot(inda,t*ones(1,length(inda)),'r.')
    end
    
    tline = fgetl(fid);
end

figure(2)
fname = sprintf('%s_tcells.pdf',fname_base);
print(fname,'-dpdf')


fclose(fid);



