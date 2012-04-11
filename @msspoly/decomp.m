function [x,p,M]=decomp(q)
% function [x,p,M]=decomp(q)
%
% If q=const then x is empty 0-by-1 msspoly, p=zeros(0,1), M=double(p(:)).
% Otherwise x is k-by-1 free msspoly, p is a sparse N-by-k nonnegative 
% integer (all rows of p are different), M is (m*n)-by-N sparse, such that
%   double(subs(q(:),x,z))=M*prod(repmat(z',N,1).^p,2),size(q)).

[ms,ns]=size(q.s);
if ns==3,
    x=msspoly(zeros(0,1));
    p=zeros(0,1);
    M=sparse(q.s(:,1),q.s(:,2),q.s(:,3));
    return
end
k=(ns-3)/2;
id=q.s(:,3:2+k);
dd=q.s(:,3+k:2+2*k);
cc=q.s(:,ns);

idx=mss_unique(id(:));     % list of unique variable id's, including 0
idx0=idx(idx~=0);          % list of unique variable id's, excluding 0
nx=length(idx0);
x=msspoly(nx,1,[(1:nx)' ones(nx,1) idx0 ones(nx,2)]);

idm=mss_match(idx0,id);
ii=repmat((1:ms)',1,k);
ee=(idm>0);
idm=idm(ee);
ddm=dd(ee);
iim=ii(ee);
p=zeros(ms,nx);
p(iim+ms*(idm-1))=ddm;
M=sparse(q.s(:,1)+q.m*(q.s(:,2)-1),(1:ms)',q.s(:,ns),q.m*q.n,ms);
[p,s]=mss_unique(p);
M=M*s;