function [A,B,C,K,S]=mss_nofree(a,b,c,k)
% function [A,B,C,K,S]=mss_nofree(a,b,c,k)
%
% remove free variables from a SeDuMi program a*x=b, c*x->min to
% non-negative ones in A*X=B, C*X->min, where x=S*X

[m,n]=size(a);
A=a; B=b; C=c; K=k; S=speye(n);
if ~isfield(K,'f'), return; end
r=K.f;
if r==0, return; end
A=[A(:,1:r) -A(:,1:r) A(:,r+1:n)];
C=[C(1:r);-C(1:r);C(r+1:n)];
S=[speye(r) -speye(r) sparse([],[],[],r,n-r); ...
    sparse([],[],[],n-r,2*r) speye(n-r)];
K=rmfield(K,'f');
if isfield(K,'l'),
    K.l=K.l+2*r;
else
    K.l=2*r;
end