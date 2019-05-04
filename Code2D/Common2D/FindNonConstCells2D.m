function nconst_ind = FindNonConstCells2D(Q)

F_EPS      = 1e-2;
Mval       = max(Q); 
mval       = min(Q);
nconst_ind = find(abs(Mval-mval) > F_EPS*max(abs(Mval),abs(mval)));

return