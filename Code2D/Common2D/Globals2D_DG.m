% Purpose: declare global variables

% global Np Nfp N K
% global r s 
% global Dr Ds LIFT Drw Dsw MassMatrix
% global Fx Fy nx ny jac Fscale J
% global vmapM vmapP vmapB mapB Fmask
% global BCTag mapBC_list vmapBC_list GEBC_list
% global rx ry sx sy J sJ
% global rk4a rk4b rk4c
% global Nfaces EToE EToF EToV
% global BFaces PerBToB_map PerBFToF_map UseMeshPerData PShift
% global V invV
% global x y NODETOL VX VY
% global FinalTime model 
% global tstamps % excluding the initial time. Must be >= 1.
% global AVG2D AVG1D_1 AVG1D_2 AVG1D_3 % Averaging matrices
% global ProjectFromNb2D CK % Needed for Fu-Shu indicator
% global Indicator Limiter ind_var lim_var
% global patch_alphas TVBM TVBnu dx CFL fixed_dt dx2
% global EToGE KG xG yG MMAP BC_flags BC_ess_flags 
% global facemid1  facemid2  facemid3
% 
% In = 1; Out = 2; Slip = 3; Far = 4; Dirichlet = 5; Sym = 6; Periodic=7;

% Low storage Runge-Kutta coefficients
% rk4a = [            0.0 ...
%         -567301805773.0/1357537059087.0 ...
%         -2404267990393.0/2016746695238.0 ...
%         -3550918686646.0/2091501179385.0  ...
%         -1275806237668.0/842570457699.0];
% rk4b = [ 1432997174477.0/9575080441755.0 ...
%          5161836677717.0/13612068292357.0 ...
%          1720146321549.0/2090206949498.0  ...
%          3134564353537.0/4481467310338.0  ...
%          2277821191437.0/14882151754819.0];
% rk4c = [             0.0  ...
%          1432997174477.0/9575080441755.0 ...
%          2526269341429.0/6820363962896.0 ...
%          2006345519317.0/3224310063776.0 ...
%          2802321613138.0/2924317926251.0]; 
