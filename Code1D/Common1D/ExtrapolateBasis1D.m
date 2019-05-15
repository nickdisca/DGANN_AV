% function [extrap_left1D, extrap_right1D] = ExtrapolateBasis1D(N,r,shift)
% 
% % Purpose : Get extrapolated modal basis values on once cell to the left 
% % and one to the right 
% 
% extrap_left1D  = zeros(length(r),N+1);
% extrap_right1D = zeros(length(r),N+1);
% for j=1:N+1
%     extrap_left1D(:,j)  = JacobiP(r(:)-shift, 0, 0, j-1);
%     extrap_right1D(:,j) = JacobiP(r(:)+shift, 0, 0, j-1);
% end
% return
