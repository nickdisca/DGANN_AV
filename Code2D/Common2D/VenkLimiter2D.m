function Qlim = VenkLimiter2D(Q)         
                                         
Globals2D_DG;

% NOT VERY HAPPY WITH THIS LIMITER YET!!!!!

% Extracting linear part of solution
% We only keep the modes 1,2 and N+2
if(N==1)
    Ql = Q;
else
    Qm              = invV*Q;
    killmodes       = [3:N+1,N+3:Np];
    Qm(killmodes,:) = 0;
    Ql              = V*Qm;
end

Qlim = Ql;

% Find neighbors in patch
E1 = EToE(:,1)'; E2 = EToE(:,2)'; E3 = EToE(:,3)';

% Get cell averages of patch
Qavg0  = AVG2D*Ql; Qavg1 = Qavg0(E1); Qavg2 = Qavg0(E2); Qavg3 = Qavg0(E3); 
MaxAvg = max([Qavg0;Qavg1;Qavg2;Qavg3],[],1); 
MinAvg = min([Qavg0;Qavg1;Qavg2;Qavg3],[],1); 


% Finding alphas of Venkatakrishnan limiter
phi        = @(y) (y.^2 + 2*y)./(y.^2 + y + 2);
AvgArray   = ones(3*Nfp,1)*Qavg0;
MaxArray   = ones(3*Nfp,1)*MaxAvg;
MinArray   = ones(3*Nfp,1)*MinAvg;
Den        = Q(Fmask(:),:) - AvgArray;
ind1 = find(Den > 0);  
ind2 = find(Den < 0);  
ind3 = find(Den == 0); 
Alp1 = 0*Den; Alp2 = 0*Den; Alp3 = 0*Den;

Alp1(ind1) = phi((MaxArray(ind1) - AvgArray(ind1))./Den(ind1));
Alp2(ind2) = phi((MinArray(ind2) - AvgArray(ind2))./Den(ind2));
Alp3(ind3) = 1;

AlphaArray = Alp1 + Alp2 + Alp3;
AlphaMin   = min(AlphaArray,[],1);


% Only limit gradients if AplhaMin<1.
indl            = find(AlphaMin < 1);

% Find Gradient of cells with AlphaMin < 1
Qr = Dr*Ql(:,indl); 
Qs = Ds*Ql(:,indl);
Qx = rx(:,indl).*Qr + sx(:,indl).*Qs;
Qy = ry(:,indl).*Qr + sy(:,indl).*Qs;

% Limit gradient
Qlim(:,indl) = ones(Np,1)*Qavg0(indl) ...
               + ((ones(Np,1)*AlphaMin(indl))).*(((eye(Np) - AVG2D)*x(:,indl)).*Qx...
                                             + ((eye(Np) - AVG2D)*y(:,indl)).*Qy);                                         

return;


