function Qlim = MRSLimiter2D(Q,QG,ind)

Globals2D_DG;

alphah =@(h) 0*h.^(1.5);
phif   =@(y) min(abs(y)/1.1,ones(size(y)));

Qlim = Q;
if(~isempty(ind))
       
    % Find neighbors in patch
    E1 = EToE(:,1)'; E2 = EToE(:,2)'; E3 = EToE(:,3)';
    
    % Get cell averages
    Qavg0  = AVG2D*Q; 
    
    % Get max and minimum values along traingle edges
    Qmax0  = max(Q);  Qmin0 = min(Q);
    Qmax1  = Qmax0(E1); Qmax2 = Qmax0(E2); Qmax3 = Qmax0(E3);
    Qmin1  = Qmin0(E1); Qmin2 = Qmin0(E2); Qmin3 = Qmin0(E3);
    
    
    % Replacing boundary element neighbours with ghost neighbours
    QGmax = max(QG);  QGmin = min(QG);
    GE1 = find(EToGE(:,1))';  GE2 = find(EToGE(:,2))'; GE3 = find(EToGE(:,3))';
    Qmax1(GE1) = QGmax(EToGE(GE1,1)); Qmax2(GE2) = QGmax(EToGE(GE2,2)); Qmax3(GE3) = QGmax(EToGE(GE3,3));
    Qmin1(GE1) = QGmin(EToGE(GE1,1)); Qmin2(GE2) = QGmin(EToGE(GE2,2)); Qmin3(GE3) = QGmin(EToGE(GE3,3));

    % Switching to only flaggd cells
    Qavg0  = Qavg0(ind);
    Qmax0  = Qmax0(ind); Qmax1  = Qmax1(ind); Qmax2  = Qmax2(ind); Qmax3  = Qmax3(ind);
    Qmin0  = Qmin0(ind); Qmin1  = Qmin1(ind); Qmin2  = Qmin2(ind); Qmin3  = Qmin3(ind);
    dxl    = dx(ind);
    alphaval = alphah(dxl);
    
    % Finding MaxVal and MinVal
    MaxVal = max([Qavg0 + alphaval;Qmax1;Qmax2;Qmax3],[],1);
    MinVal = min([Qavg0 - alphaval;Qmin1;Qmin2;Qmin3],[],1);
    
    % Finding Thetamax and Thetamin
    eps = 1.0e-12;
    Thetamax = phif((MaxVal-Qavg0)./(Qmax0 - Qavg0));
    Thetamin = phif((MinVal-Qavg0)./(Qmin0 - Qavg0));
    Thetaval = min([ones(size(MaxVal));Thetamax;Thetamin],[],1);
    
    Qlim(:,ind) = ones(Np,1)*(Qavg0.*(1-Thetaval)) + (ones(Np,1)*Thetaval).*Qlim(:,ind);
    
end

return;
