function P=sparse(p)
% function P=sparse(p)
%
% for an msspoly matrix p=p(y) produces sparse matrix P of
% minimal dimension such that p(y)x=0 for all y iff Px=0

% AM 11.01.09

if isempty(p.s), P=sparse([],[],[],0,p.n); return; end
[ms,ns]=size(p.s);
if ms==1, P=sparse(1,p.s(1,2),p.s(1,ns),1,p.n); return; end
[ss,ii]=sortrows([p.s(:,1) p.s(:,3:ns-1)]);
jj=p.s(ii,2);
cc=p.s(ii,ns);
ee=cumsum([1;1-all(ss(1:ms-1,:)==ss(2:ms,:),2)]);
P=sparse(ee,jj,cc,ee(length(ee)),p.n);