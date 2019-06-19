function [mu_piece] = Euler1D_viscosity(q,qold,Viscosity,Problem,Mesh,dt,iter, Net)
% Return the elementwise values for the viscosity - the
% smoothing is done afterwards

%compute maximum wave speed 
pre    = (Problem.gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
c_sound = sqrt(abs(Problem.gas_gamma*pre./q(:,:,1)));
lambda = c_sound + abs(q(:,:,2)./q(:,:,1));
local_wave_sp=max(lambda);

%select viscosity variable
if strcmp(Viscosity.visc_var,'density')
    var=q(:,:,1);
else
    error('Viscosity has to be computed using density');
end

%compute viscosity according to the selected model
switch Viscosity.model
    
    case 'NONE'
        mu_piece=zeros(1,Mesh.K);
        
        
    case "EV"
        [entropy,ent_flux] = Set_Euler_entropy1D(Problem.gas_gamma);
        
        E_old=entropy(qold); E_new=entropy(q);
        F_old=ent_flux(qold); F_new=ent_flux(q);
        % f_prime_old=dflux(uold); f_prime_new=dflux(u);
        
        mu_max=Viscosity.c_max*Mesh.hK/Mesh.N.*local_wave_sp;
        
        R=(E_new-E_old)/dt+(Mesh.Dr*F_old+Mesh.Dr*F_new)./Mesh.hK;
        if (iter==0)
            R=zeros(size(R));
        end
        
        % E_ext = Apply_BC1D(E_new,{'N',0.0,'N',0.0});
        F_ext = Apply_BC1D(F_new(Mesh.VtoE),{'N',0.0,'N',0.0});
        % f_prime_ext=extendDG(f_prime_new,{'N',0.0,'N',0.0});
        cL=F_ext(2,1:Mesh.K)-F_ext(1,2:Mesh.K+1); 
        cR=F_ext(2,2:Mesh.K+1)-F_ext(1,3:Mesh.K+2);
        J=max( abs(cL),abs(cR) ) ./ ( Mesh.hK/Mesh.N );
        
        Norm_E=1;
        
        DD=max(max(abs(R)),abs(J))./Norm_E;
        mu_E=Viscosity.c_E.*(Mesh.hK/Mesh.N).^2.*DD; 
        mu_piece=min(mu_E,mu_max);   
        
        
    case "MDH"
        mu_max=Viscosity.c_max*Mesh.hK/Mesh.N.*local_wave_sp;
        var_modal=Mesh.invV*var;
        S_K=var_modal(Mesh.N+1,:).^2./sum(var_modal.^2);
        s_K=log10(S_K);
        s0=-(Viscosity.c_A+4*log10(Mesh.N));
        mu1=((s0-Viscosity.c_k)<=s_K).*(s_K<=(s0+Viscosity.c_k)); mu2=s_K>(s0+Viscosity.c_k);
        mu_piece=mu_max.*(mu1/2.*(1+sin(pi/(2*Viscosity.c_k)*(s_K-s0))) + mu2);
        mu_piece(isnan(mu_piece))=0;
 
        
    case "MDA"
        mu_max=Viscosity.c_max*Mesh.hK/Mesh.N.*local_wave_sp;
        
        if (Mesh.N<=2)
            mu_piece=mu_max;
        else
        
            var_modal=Mesh.invV*var;
        
            b_i=repmat([1:Mesh.N]'.^(-Mesh.N)./sqrt(sum([1:Mesh.N]'.^(-2*Mesh.N))),1,Mesh.K);
        
            norm2=repmat(sum(var_modal(:,:).^2.*Mesh.hK/2),Mesh.N,1);
            var_breve=sqrt(var_modal(2:end,:).^2+norm2.*(b_i.^2));
        
            var_check=zeros(Mesh.N,Mesh.K);
            for ii=1:Mesh.N
                var_check(ii,:)=max(abs(var_breve(min(ii,Mesh.N-1):end,:)),[],1);
            end
        
            A11=Mesh.N; A12=-sum(log(1:Mesh.N)); A21=A12; A22=sum(log(1:Mesh.N).^2);
            B1=sum(log(abs(var_check(:,:)))); B2=-sum(log(abs(var_check(:,:))).*repmat(log(1:Mesh.N)',1,Mesh.K));
            coef=[A11 A12; A21 A22]\[B1; B2];
            coef(isnan(coef))=0;
                        
            mu1=(coef(2,:)<1); mu2=(coef(2,:)>=1).*(coef(2,:)<=3);
            mu_piece=mu_max.*(mu1 + mu2.*(1-(coef(2,:)-1)/2));
        
        end           

    
    case "NN"
        
        mu_piece = visc_MLP1D(var,Net,var,Mesh.hK,local_wave_sp,Mesh.VtoE, Problem.bc_cond);
        
        c=1.0; mu_piece = c*mu_piece;
        
        
    otherwise
        error('Unknown viscosity model');
end
    
end

