function Qlim = BarthJespLimiter2D(Q,ind)

Globals2D_DG;

if(isempty(ind))
    Qlim = Q;
else
    % Extracting linear part of solution for elements marked using indicator
    % These are listed in "ind"
    % We only keep the modes 1,2 and N+2
    if(N==1)
        Ql = Q;
    else
        Qm              = invV*Q;
        killmodes       = [3:N+1,N+3:Np];
        Qm(killmodes,ind) = 0;
        Ql              = V*Qm;
    end
    Qlim = Ql;
    
    
    % Working with only the flagged elements
    % Find neighbors in patch
    E1 = EToE(ind,1)'; E2 = EToE(ind,2)'; E3 = EToE(ind,3)';
    
    % Get cell averages of patch
    Qavg0  = AVG2D*Ql; Qavg1 = Qavg0(E1); Qavg2 = Qavg0(E2); Qavg3 = Qavg0(E3);
    Qavg0  = Qavg0(ind);
    Ql     = Ql(:,ind);
    MaxAvg = max([Qavg0;Qavg1;Qavg2;Qavg3],[],1);
    MinAvg = min([Qavg0;Qavg1;Qavg2;Qavg3],[],1);
    
    % Get local x,y node coordinates, and local jacobian factors
    xl  = x(:,ind);  yl  = y(:,ind);
    rxl = rx(:,ind); sxl = sx(:,ind);  ryl = ry(:,ind); syl = sy(:,ind);
    
    
    
    % Finding alphas of Barth Jesperson limiter
    AvgArray   = ones(3*Nfp,1)*Qavg0;
    MaxArray   = ones(3*Nfp,1)*MaxAvg;
    MinArray   = ones(3*Nfp,1)*MinAvg;
    Num1       = MaxArray - AvgArray;
    Num2       = MinArray - AvgArray;
    Den        = Ql(Fmask(:),:) - AvgArray;
    
    ind1 = find(Ql(Fmask(:),:) - MaxArray > 0);
    ind2 = find(Ql(Fmask(:),:) - MinArray < 0);
    
    AlphaArray       = ones(size(Num1));
    AlphaArray(ind1) = Num1(ind1)./Den(ind1);
    AlphaArray(ind2) = Num2(ind2)./Den(ind2);
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
    assert(min(min(Ql - ones(Np,1)*MinAvg))>-eps)
    assert(min(min(-Ql + ones(Np,1)*MaxAvg))>-eps)
    % Updating Qlim
    Qlim(:,ind) = Ql;
    
end

return;
