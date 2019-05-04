function ind = KXRCF2D_ncon(Q,QG)

% Purpose: find all the troubled-cells for variable Q using KXRCF 
% indicator. We exclude almost constant cells.

Globals2D_DG;

eps = 1e-10;

% Finding non-constant cells
F_EPS    = 1e-2;
Mval     = max(Q); 
mval     = min(Q);
ncon_ind = find(abs(Mval-mval) > F_EPS*max(abs(Mval),abs(mval)));

% Computing '-' and '+' traces
QM(:,:) = Q(vmapM); QP(:,:) = Q(vmapP);

% set boundary conditions by modifying positive traces
if(~isempty(mapBC_list))
    QP(mapB) = QG(Fmask(:,Nfaces),:);
end



% Find neighbors in patch
E1 = EToE(:,1)'; E2 = EToE(:,2)'; E3 = EToE(:,3)';

% Extracting linear part of solution for elements. 
% We only keep the modes 1,2 and N+2
if(N==1)
    Ql  = Q;
    QGl = QG;
else
    Qm              = invV*Q;
    killmodes       = [3:N+1,N+3:Np];
    Qm(killmodes,:) = 0;
    Ql              = V*Qm;
    
    QGm              = invV*QG;
    killmodes        = [3:N+1,N+3:Np];
    QGm(killmodes,:) = 0;
    QGl              = V*QGm;
end

% Get cell averages of patch. For the time being, Neuman BC used (except 
% for periodic BC)
AVG0 = AVG2D*Ql;
AVGn = [AVG0(E1); AVG0(E2); AVG0(E3)];

% Replacing boundary element neighbours with ghost neighbours
AVGGE = AVG2D*QGl;
GE1 = find(EToGE(:,1))';
GE2 = find(EToGE(:,2))';
GE3 = find(EToGE(:,3))';
AVGn(1,GE1) = AVGGE(EToGE(GE1,1)); 
AVGn(2,GE2) = AVGGE(EToGE(GE2,2)); 
AVGn(3,GE3) = AVGGE(EToGE(GE3,3));

% xavg = AVG2D*x;
% yavg = AVG2D*y;
% 
% figure(100)
% subplot(2,2,1)
% scatter3(xavg,yavg,AVG0,'o')
% subplot(2,2,2)
% scatter3(xavg,yavg,AVG0,'o')
% subplot(2,2,3)
% scatter3(xavg,yavg,AVG0,'o')
% subplot(2,2,4)
% scatter3(xavg,yavg,AVG0,'o')



% Get value at face mid point value for linear projection
QCF1 = facemid1*Ql(Fmask(:,1),:); 
QCF2 = facemid2*Ql(Fmask(:,2),:); 
QCF3 = facemid3*Ql(Fmask(:,3),:);

% Qtildes
Qtilde1 = QCF1 - AVG0; Qtilde2 = QCF2 - AVG0; Qtilde3 = QCF3 - AVG0;

% DelQs
ind1    = reshape(patch_alphas(1,3,:),[1,K]); ind1xy  = sub2ind(size(AVGn),ind1,1:K);
ind2    = reshape(patch_alphas(2,3,:),[1,K]); ind2xy  = sub2ind(size(AVGn),ind2,1:K);
ind3    = reshape(patch_alphas(3,3,:),[1,K]); ind3xy  = sub2ind(size(AVGn),ind3,1:K);
DelQ1   = sum( reshape(patch_alphas(1,1:2,:),[2,K]).* ...
               [AVGn(1,:) - AVG0; AVGn(ind1xy) - AVG0] );
DelQ2   = sum( reshape(patch_alphas(2,1:2,:),[2,K]).* ...
               [AVGn(2,:) - AVG0; AVGn(ind2xy) - AVG0] );
DelQ3   = sum( reshape(patch_alphas(3,1:2,:),[2,K]).* ...
               [AVGn(3,:) - AVG0; AVGn(ind3xy) - AVG0] );           

% Get limited slopes/modes           
s1 = minmodB([Qtilde1;TVBnu*DelQ1],TVBM,dx);
s2 = minmodB([Qtilde2;TVBnu*DelQ2],TVBM,dx);
s3 = minmodB([Qtilde3;TVBnu*DelQ3],TVBM,dx);

% Fix slope if s1+s2+s3 = 0
ind0   = find(abs(s1+s2+s3)>eps);
% length(ind0)
pos    = max(0,s1(ind0)) + max(0,s2(ind0)) + max(0,s3(ind0));
neg    = max(0,-s1(ind0)) + max(0,-s2(ind0)) + max(0,-s3(ind0));
thetap = min(1,neg./pos); 
thetam = min(1,pos./neg);
s1_fix = thetap.*max(0,s1(ind0)) - thetam.*max(0,-s1(ind0));
s2_fix = thetap.*max(0,s2(ind0)) - thetam.*max(0,-s2(ind0));
s3_fix = thetap.*max(0,s3(ind0)) - thetam.*max(0,-s3(ind0));
% max(abs(s1_fix - s1(ind0)))
% max(abs(s2_fix - s2(ind0)))
% max(abs(s3_fix - s3(ind0)))
s1(ind0) = s1_fix;
s2(ind0) = s2_fix;
s3(ind0) = s3_fix;


% Get limited face mid-point values
QCFL1 = AVG0 + s1; QCFL2 = AVG0 + s2; QCFL3 = AVG0 + s3;

% Mark cells where even one of the face values has changed 
% QCF1_true = facemid1*Q(Fmask(:,1),:);
% QCF2_true = facemid2*Q(Fmask(:,2),:); 
% QCF3_true = facemid3*Q(Fmask(:,3),:);


% for i = 1:K
%     figure(1001)
%     plot(0,QCF1_true(i),'o',0,QCF1(i),'x',0,QCFL1(i),'s')
%     
%     figure(1002)
%     plot(0,QCF2_true(i),'o',0,QCF2(i),'x',0,QCFL2(i),'s')
%     
%     figure(1003)
%     plot(0,QCF3_true(i),'o',0,QCF3(i),'x',0,QCFL3(i),'s')
%     
%     pause(0.1)
% end

%ind = find(abs(QCFL1 - QCF1_true) > eps | abs(QCFL2 - QCF2_true) > eps | abs(QCFL3 - QCF3_true) > eps);
ind = find(abs(QCFL1 - QCF1) > eps | abs(QCFL2 - QCF2) > eps | abs(QCFL3 - QCF3) > eps);


% Only keeping non-constant cells
ind = intersect(ind,ncon_ind);



return
