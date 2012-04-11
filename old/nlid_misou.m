function [f,h,r,U,Q,D]=nlid_misou(Z,p)
% function [f,h,r,U,Q,D]=nlid_misou(Z,p)
%
% Nonlinear System Identification for 
% Multiple Input Single Output systems: the unit cube case
%
% Generate incrementally stable polynomial MISO model 
% (where input w and output v have values in the interval (-1,1)):
%       f(v(t),v(t-1),...,v(t-n),w0(t),...,wk(t))=0           (*)
% from data 
%       Z(i,:)=[v(t(i)) v(t(i)-1) ... v(t(i)-n) w0(t(i)) ... wk(t(i))]
%
% INPUTS:
%   Z  -   mz-by-(n+k+2) real: i/o data, all entries within the unit cube
%   p  -   mp-by-(k+2) integer>0: determines which monomials to use
% 
% OUTPUTS:
%   f   -  scalar msspoly in z=[v;w]: v=[v0;...;vn], w=[w0,...,wk]
%   h   -  scalar msspoly in v(2:n+1), d(2:n+1), where d=[d0;...;dn]
%   r   -  scalar msspoly in z
%   U   -  column of msspoly monomials in variables z,d
%   Q   -  symmetric positive semidefinite matrix
% 
% h,r,U,Q certify incremental stability of (*) via the polynomial identity
%   r+d(1)*(f+diff(f,v,d))-d(1)^2+h-subs(h,x0,x1)=U'*Q*U
% where x0=[v(2:n+1);d(2:n+1)], x1=[v(1:n);d(1:n)].
% Robust fidelity measure sum_i r(Z(i,:)) is minimized.
%
% p  defines U=[U0;U1*d(1);U2*d(2);...;Un*d(n)], where 
% U0 has monomials (v.^a)*(w.^b): sum(a)/p(i,1)+sum(b./p(i,2:k+2))<=1
% U1 has monomials (v.^a)*(w.^b): (sum(a)+1)/p(i,1)+sum(b./p(i,2:k+2))<=1
% Uk (k>1) has monomials v(k:n+1).^a: sum(a)<max(p(:,1))

if nargin<2, error('2 inputs required'); end
if ~isa(Z,'double'), error('input 1 not a "double"'); end
if any(abs(Z(:))>1), error('input 1 out of bounds'); end
[mz,nz]=size(Z);
if mz*nz==0, error('input 1 is empty'); end
if ~mint_isint(p), error('input 2 not integer'); end
if any(p(:)<1), error('input 2 not positive'); end
[mp,np]=size(p);
if mp*np==0, error('input 2 is empty'); end
k=np-2;
if k<0, error('input 2 too short'); end
n=nz-k-2;
if n<1, error('input 1: not enough columns'); end

v=msspoly('v',n+1);  
w=msspoly('w',k+1); 
d=msspoly('d',n+1);
z=[v;w];

K0=zeros(1+nz*mp,nz);
for i=1:mp, 
    K0(1+nz*(i-1):nz+i,:)=diag(p(i,:)); 
end
U0=recomp(z,mint_ch(K0));



if r>1,
    vd0=[v(2:r);d(2:r)];
    vd1=[v(1:r-1);d(1:r-1)];
    if ~isfunction(Z,vd0), error('illegal variables in the 6th input'); end
end
fprintf('\nnlid_miso: pre-computing.');tic;
C0=zeros(mr,1); 
for i=1:T, 
    C0=C0+double(subs(R,vw,[vv(i,:) ww(i,:)]')); 
    if mod(i,100)==0, fprintf('.'); end
end
toc
pr=mssprog;
[pr,a]=new(pr,mf,'free','a');
[pr,b]=new(pr,mr,'free','b');
[pr,Q]=new(pr,mu,'psd','Q');
f=a'*F;
g=f+diff(f,v,d);
e=b'*R+2*d(1)*g-d(1)^2-U'*Q*U;
if r>1,
    [pr,P]=new(pr,mz,'psd','P');
    h0=Z'*P*Z ;
    h1=subs(h0,vd0,vd1);
    e=e+h0-h1;
end
pr=eq(pr,e);
pr=sedumi(pr,C0'*b);
a=double(pr,a);
b=double(pr,b);
Q=double(pr,Q);
if r>1,
    P=double(pr,P);
else
    P=1;
end
fprintf('nlid_miso: upper bound %e, certified: %f>0, %f>0\n', ...
    C0'*b,min(eig(Q)),min(eig(P)))

