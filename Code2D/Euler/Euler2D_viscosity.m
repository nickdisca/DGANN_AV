function [mu_piece] = Euler2D_viscosity(Q,Qold, local_wave_sp,Viscosity,Problem,Mesh,dt,iter, Net)
% Return the elementwise values for the viscosity - the
% smoothing is done afterwards

%compute maximum wave speed 
local_wave_sp=max(local_wave_sp);

%select viscosity variable
if strcmp(Viscosity.visc_var,'density')
    var=Q(:,:,1);
else
    error('Viscosity has to be computed using density');
end

%compute viscosity according to the selected model
switch Viscosity.model
    
    case 'NONE'
        mu_piece=zeros(1,Mesh.K);
        
        
    case "EV"

        [E_old,F_old_X,F_old_Y] = Compute_Euler_entropy2D(Qold,Problem.gas_gamma);
        [E_new,F_new_X,F_new_Y] = Compute_Euler_entropy2D(Q,Problem.gas_gamma);
        
        mu_max=Viscosity.c_max*2*Mesh.dx2/Mesh.N.*local_wave_sp;
        
        divFold = Mesh.rx.*(Mesh.Dr*F_old_X) + Mesh.sx.*(Mesh.Ds*F_old_X) + Mesh.ry.*(Mesh.Dr*F_old_Y) + Mesh.sy.*(Mesh.Ds*F_old_Y);
        divFnew = Mesh.rx.*(Mesh.Dr*F_new_X) + Mesh.sx.*(Mesh.Ds*F_new_X) + Mesh.ry.*(Mesh.Dr*F_new_Y) + Mesh.sy.*(Mesh.Ds*F_new_Y);
        
        R=(E_new-E_old)/dt+(divFold+divFnew)/2;
        if (iter==0)
            R=zeros(size(R));
        end
        
        
        FXP = F_new_X(Mesh.vmapP); FXM=F_new_X(Mesh.vmapM); FYP = F_new_Y(Mesh.vmapP); FYM=F_new_Y(Mesh.vmapM); 
        
        Jump=max(abs(reshape(FXP-FXM,(Mesh.N+1)*3,Mesh.K).*Mesh.nx+reshape(FYP-FYM,(Mesh.N+1)*3,Mesh.K).*Mesh.ny)./(ones((Mesh.N+1)*3,1)*2*Mesh.dx2/Mesh.N));
        
        Norm_E=1;
        
        DD=max(max(abs(R)),abs(Jump))/Norm_E;
        mu_E=Viscosity.c_E*(2*Mesh.dx2/Mesh.N).^2.*DD; 
        mu_piece=min(mu_E,mu_max);  

        
        
    case "MDH"
        
        mu_max=Viscosity.c_max*2*Mesh.dx2/Mesh.N.*local_wave_sp;
        
        var_modal=Mesh.invV*var;
        ids_high=(1:(Mesh.N+1))*(Mesh.N+1)-(1:Mesh.N+1).*(0:Mesh.N)/2; 
        S_K=sum(var_modal(ids_high,:).^2)./sum(var_modal.^2);
        s_K=log10(S_K);
        s0=-(Viscosity.c_A+4*log10(Mesh.N));
        mu1=((s0-Viscosity.c_k)<=s_K).*(s_K<=(s0+Viscosity.c_k)); mu2=s_K>(s0+Viscosity.c_k);
        mu_piece=mu_max.*(mu1/2.*(1+sin(pi/(2*Viscosity.c_k)*(s_K-s0))) + mu2);
        mu_piece(isnan(mu_piece))=0;

 
        
    case "MDA"
        
        mu_max=Viscosity.c_max*2*Mesh.dx2/Mesh.N.*local_wave_sp;
        
        if (Mesh.N<=2)
            mu_piece=mu_max;
        else
            
            var_faces=reshape(var(Mesh.vmapM),3*(Mesh.N+1),Mesh.K);
            Vand1D=Vandermonde1D(Mesh.N,JacobiGL(0,0,Mesh.N));
            
            for face=1:Mesh.Nfaces
                
                ids_face=(Mesh.N+1)*(face-1)+1:(Mesh.N+1)*face;
                var_face=var_faces(ids_face,:);
                var_face_modal=Vand1D\var_face;
                
                b_i=repmat([1:Mesh.N]'.^(-Mesh.N)./sqrt(sum([1:Mesh.N]'.^(-2*Mesh.N))),1,Mesh.K);
                
                J_bd=Mesh.Fscale(ids_face,:).*(ones(Mesh.N+1,1)*Mesh.J(1,:));
                norm2=repmat(sum(var_face_modal.^2).*J_bd(1,:),Mesh.N,1);
                var_face_breve=sqrt(var_face_modal(2:end,:).^2+norm2.*(b_i.^2));
                
                for ii=1:Mesh.N
                    var_face_check(ii,:)=max(abs(var_face_breve(min(ii,Mesh.N-1):end,:)),[],1);
                end
                
                A11=Mesh.N; A12=-sum(log(1:Mesh.N)); A21=A12; A22=sum(log(1:Mesh.N).^2);
                B1=sum(log(abs(var_face_check(:,:)))); B2=-sum(log(abs(var_face_check(:,:))).*repmat(log(1:Mesh.N)',1,Mesh.K));
                coef=[A11 A12; A21 A22]\[B1; B2];
                coef(isnan(coef))=0;
                
                estimated_decay(face,:)=coef(2,:);
            end
            
            global_decay=max(estimated_decay);
            
            mu1=(global_decay<1); mu2=(global_decay>=1).*(global_decay<=3);
            mu_piece=mu_max.*(mu1 + mu2.*(1-(coef(2,:)-1)/2));
            
            
        end

       
  
    
    case "NN"
        
        mu_piece = visc_MLP2D(var,Net,var,2*Mesh.dx2,local_wave_sp, Mesh);
        
        c=1.0; mu_piece = c*mu_piece;
        
        
    otherwise
        error('Unknown viscosity model');
end
    
end