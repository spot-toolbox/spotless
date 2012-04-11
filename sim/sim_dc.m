function [h,u,f]=sim_dc(x,d)
% function [h,u,f]=sim_dc(x,d)
%
% INPUTS:
%   x  -  m-by-n real
%   d  -  1-by-n positive integer
% OUTPUT:
%   h  -  non-negative scalar
%   u  -  normalized and transformed x
%   f  -  function handle

x=real(double(x));
d=max(1,round(real(d(1,:))));
[m,n]=size(x);
if n~=length(d), error('incompatible dimensions'); end
x0=min(x,[],1);
x1=max(x,[],1);
s=(x1>x0);
if ~any(s), h=0; return; end   % a single point
x0=x0(s);
x1=x1(s);
u=acos((x(:,s)-repmat(x0,m,1))./repmat(x1-x0,m,1));
dd=mint_ch(mint_down(d));
dd=dd(2:size(dd,1),:);
M=u*dd';
e=repmat(1/sqrt(2),m,1);
M=[e cos(M) sin(M)];
M=M'*M;
if nargout<3,
    r=min(eig(M));
else
    [V,D]=eig(M);
    [r,i]=min(diag(D));
    v=V(:,i);
    f=@(w)[repmat(1/sqrt(2),size(w,1),1) cos(w*dd') sin(w*dd')]*v;
end    
h=sqrt(r/m);