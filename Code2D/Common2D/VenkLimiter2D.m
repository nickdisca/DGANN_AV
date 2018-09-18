function Qlim = VenkLimiter2D(Q,QG,ind)         

Globals2D_DG;

Qlim = Q;
if(~isempty(ind))
    % Extracting linear part of solution for elements marked using indicator
    % These are listed in "ind"
    % We only keep the modes 1,2 and N+2
    if(N==1)
        Ql = Q;
        QGl = QG;
    else
        Qm              = invV*Q;
        killmodes       = [3:N+1,N+3:Np];
        Qm(killmodes,ind) = 0;
        Ql              = V*Qm;
        
        QGm              = invV*QG;
        killmodes        = [3:N+1,N+3:Np];
        QGm(killmodes,:) = 0;
        QGl              = V*QGm;
    end
       
    % Find neighbors in patch
    E1 = EToE(:,1)'; E2 = EToE(:,2)'; E3 = EToE(:,3)';
    
    % Get cell averages of patch
    Qavg0  = AVG2D*Ql; Qavg1 = Qavg0(E1); Qavg2 = Qavg0(E2); Qavg3 = Qavg0(E3);
    
    
    % Replacing boundary element neighbours with ghost neighbours
    QGavg = AVG2D*QGl;
    GE1 = find(EToGE(:,1))';  GE2 = find(EToGE(:,2))'; GE3 = find(EToGE(:,3))';
    Qavg1(GE1) = QGavg(EToGE(GE1,1));
    Qavg2(GE2) = QGavg(EToGE(GE2,2));
    Qavg3(GE3) = QGavg(EToGE(GE3,3));


    % Switching to only flaggd cells
    Qavg0  = Qavg0(ind); Qavg1  = Qavg1(ind); Qavg2  = Qavg2(ind); Qavg3  = Qavg3(ind);
    Ql     = Ql(:,ind);
    
    % Finding avg bounds
    MaxAvg = max([Qavg0;Qavg1;Qavg2;Qavg3],[],1);
    MinAvg = min([Qavg0;Qavg1;Qavg2;Qavg3],[],1);
    
    % Get local x,y node coordinates, and local jacobian factors
    xl  = x(:,ind);  yl  = y(:,ind);
    rxl = rx(:,ind); sxl = sx(:,ind);  ryl = ry(:,ind); syl = sy(:,ind);
    
    
    
    % Finding alphas of Venkatakrishnan limiter
    phi        = @(y) (y.^2 + 2*y)./(y.^2 + y + 2);
    AvgArray   = ones(3*Nfp,1)*Qavg0;
    MaxArray   = ones(3*Nfp,1)*MaxAvg;
    MinArray   = ones(3*Nfp,1)*MinAvg;
    Num1       = MaxArray - AvgArray;
    Num2       = MinArray - AvgArray;
    Den        = Ql(Fmask(:),:) - AvgArray;
    
    ind1 = find(Ql(Fmask(:),:) - MaxArray > 0);
    ind2 = find(Ql(Fmask(:),:) - MinArray < 0);
    
    AlphaArray       = ones(size(Num1));
    AlphaArray(ind1) = phi(Num1(ind1)./Den(ind1));
    AlphaArray(ind2) = phi(Num2(ind2)./Den(ind2));
    AlphaMin         = min(AlphaArray,[],1);
    

    % Only limit gradients if AplhaMin<1.
    indl            = find(AlphaMin < 1);
    
    % Limit gradient
    if(~isempty(indl))
        
        % Find Gradient of cells with AlphaMin < 1
        Qr = Dr*Ql(:,indl);
        Qs = Ds*Ql(:,indl);
        Qx = rxl(:,indl).*Qr + sxl(:,indl).*Qs;
        Qy = ryl(:,indl).*Qr + syl(:,indl).*Qs;
        
        
        Ql(:,indl) = ones(Np,1)*Qavg0(indl) ...
            + ((ones(Np,1)*AlphaMin(indl))).*(((eye(Np) - AVG2D)*xl(:,indl)).*Qx...
            + ((eye(Np) - AVG2D)*yl(:,indl)).*Qy);
    end
%     size(MinArray)
%     size(Ql)
    eps = 1.0e-10;
%     min(min(Ql - ones(Np,1)*MinAvg))
%     min(min(-Ql + ones(Np,1)*MaxAvg))
%    assert(min(min(Ql - ones(Np,1)*MinAvg))>-eps)
%    assert(min(min(-Ql + ones(Np,1)*MaxAvg))>-eps)
    % Updating Qlim
    Qlim(:,ind) = Ql;
    
end
                             

return;


