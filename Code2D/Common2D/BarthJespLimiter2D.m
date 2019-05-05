function Qlim = BarthJespLimiter2D(Q,QG,ind,Mesh)

Qlim = Q;
if(~isempty(ind))
    % Extracting linear part of solution for elements marked using indicator
    % These are listed in "ind"
    % We only keep the modes 1,2 and N+2
    if(Mesh.N==1)
        Ql = Q;
        QGl = QG;
    else
        Qm              = Mesh.invV*Q;
        killmodes       = [3:Mesh.N+1,Mesh.N+3:Mesh.Np];
        Qm(killmodes,ind) = 0;
        Ql              = Mesh.V*Qm;
        
        QGm              = Mesh.invV*QG;
        killmodes        = [3:Mesh.N+1,Mesh.N+3:Mesh.Np];
        QGm(killmodes,:) = 0;
        QGl              = Mesh.V*QGm;
    end
     
%     [TRI,xout,yout,interp] = GenInterpolators2D(N, N, x, y, invV);
%     figure(10)
%     PlotField2D(Ql,interp,TRI,xout,yout); axis tight; drawnow;
%     hold all
    
    % Find neighbors in patch
    E1 = Mesh.EToE(:,1)'; E2 = Mesh.EToE(:,2)'; E3 = Mesh.EToE(:,3)';
    
    % Get cell averages of patch
    Qavg0  = Mesh.AVG2D*Ql; Qavg1 = Qavg0(E1); Qavg2 = Qavg0(E2); Qavg3 = Qavg0(E3);
    
    
    % Replacing boundary element neighbours with ghost neighbours
    QGavg = Mesh.AVG2D*QGl;
    GE1 = find(Mesh.EToGE(:,1))';  GE2 = find(Mesh.EToGE(:,2))'; GE3 = find(Mesh.EToGE(:,3))';
    Qavg1(GE1) = QGavg(Mesh.EToGE(GE1,1));
    Qavg2(GE2) = QGavg(Mesh.EToGE(GE2,2));
    Qavg3(GE3) = QGavg(Mesh.EToGE(GE3,3));


    % Switching to only flaggd cells
    Qavg0  = Qavg0(ind); Qavg1  = Qavg1(ind); Qavg2  = Qavg2(ind); Qavg3  = Qavg3(ind);
    Ql     = Ql(:,ind);
    
    % Finding avg bounds
    MaxAvg = max([Qavg0;Qavg1;Qavg2;Qavg3],[],1);
    MinAvg = min([Qavg0;Qavg1;Qavg2;Qavg3],[],1);
    
    % Get local x,y node coordinates, and local jacobian factors
    xl  = Mesh.x(:,ind);  yl  = Mesh.y(:,ind);
    rxl = Mesh.rx(:,ind); sxl = Mesh.sx(:,ind);  ryl = Mesh.ry(:,ind); syl = Mesh.sy(:,ind);
    
    
    
    % Finding alphas of Barth Jesperson limiter
    AvgArray   = ones(3*Mesh.Nfp,1)*Qavg0;
    MaxArray   = ones(3*Mesh.Nfp,1)*MaxAvg;
    MinArray   = ones(3*Mesh.Nfp,1)*MinAvg;
    Num1       = MaxArray - AvgArray;
    Num2       = MinArray - AvgArray;
    Den        = Ql(Mesh.Fmask(:),:) - AvgArray;
    
    eps = 1.0e-12;
    ind1 = find(Ql(Mesh.Fmask(:),:) - MaxArray > eps);
    ind2 = find(Ql(Mesh.Fmask(:),:) - MinArray < -eps);
    
    AlphaArray       = ones(size(Num1));
    AlphaArray(ind1) = Num1(ind1)./Den(ind1);
    AlphaArray(ind2) = Num2(ind2)./Den(ind2);
    AlphaMin         = min(AlphaArray,[],1);
    

    % Only limit if AplhaMin<1.

    %indl            = find(AlphaMin < 1);
    indnl           = find(abs(AlphaMin-1)<eps);
    indl            = setdiff(1:length(AlphaMin),indnl);
    ind             = ind(indl);
    
    % Limit gradient
    if(~isempty(indl))
        
        % Find Gradient of cells with AlphaMin < 1
        Qr = Mesh.Dr*Ql(:,indl);
        Qs = Mesh.Ds*Ql(:,indl);
        Qx = rxl(:,indl).*Qr + sxl(:,indl).*Qs;
        Qy = ryl(:,indl).*Qr + syl(:,indl).*Qs;
        
%         xf = x(:,ind);
%         yf = y(:,ind);
%         QQ = Ql(:,indl);
%         figure(400)
%         scatter3(xf(:),yf(:),QQ(:))
%         hold all
        
        Ql(:,indl) = ones(Mesh.Np,1)*Qavg0(indl) ...
            + ((ones(Mesh.Np,1)*AlphaMin(indl))).*(((eye(Mesh.Np) - Mesh.AVG2D)*xl(:,indl)).*Qx...
            + ((eye(Mesh.Np) - Mesh.AVG2D)*yl(:,indl)).*Qy);
        
%         QQ = Ql(:,indl);
%         scatter3(xf(:),yf(:),QQ(:))
%         hold off
        %assert(min(min(Ql(:,indl) - ones(Np,1)*MinAvg(indl)))>-eps)
        %assert(min(min(-Ql(:,indl) + ones(Np,1)*MaxAvg(indl)))>-eps)
    end
%     size(MinArray)
%     size(Ql)
    
%     min(min(Ql - ones(Np,1)*MinAvg))
%     min(min(-Ql + ones(Np,1)*MaxAvg))
%    assert(min(min(Ql - ones(Np,1)*MinAvg))>-eps)
%    assert(min(min(-Ql + ones(Np,1)*MaxAvg))>-eps)
    % Updating Qlim
    %Qlim(:,ind) = Ql;
    Qlim(:,ind) = Ql(:,indl);
    
%     figure(11)
%     PlotField2D(Qlim,interp,TRI,xout,yout); axis tight; drawnow;
    
end

return;
