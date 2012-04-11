function [aa,bb,dd]=nlid_fit(x,y,q)
% function nlid_fit(x,y,q)
%
% simple minded rational fit

if (max(x(:))>1)||(max(y)>1)||(min(x(:))<0)||(min(y)<0),
    error('input data out of bounds');
end
[n,k]=size(x);
n0=round(2000^(1/k));
[x0,y0]=sim_select(x,y,repmat(n0,1,k));
%close(gcf);plot(xn(:,1),xn(:,2),'.');grid;pause
%close(gcf);plot(xn0(:,1),xn0(:,2),'.');grid;pause
%close(gcf);plot(yn);grid;pause
%close(gcf);plot(yn0);grid;pause
[aa,bb,dd,L]=nlid_uw2ab(y0,acos(x0),q,1);
y0h=nlid_abdx2u(aa,bb,dd,x0);
close(gcf);plot(y0,y0h,'.'); grid; pause
yh=nlid_abdx2u(aa,bb,dd,x);
close(gcf);plot(y,yh,'.'); grid; pause
fprintf('  error: bound=%f, select=%f, all=%f\n', ...
    L,sqrt(mean((y0-y0h).^2)),sqrt(mean((y-yh).^2)))