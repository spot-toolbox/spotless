function [f,r,h,U0,U1,U2,Q]=nlid_miso1(Z,p)
% function [f,r,h,U0,U1,U2,Q]=nlid_miso1(Z,p)
%
% Non-Linear System Identification for 
% Multiple Input Single Output systems (first order version)
%
% Generate implicit polynomial MISO model (input w, output v)
%       f(v(t),v(t-1),w0(t),...,wk(t))=0           (*)
% from data Z(i,:)=[v(t(i)) v(t(i)-1) w0(t(i)) ... wk(t(i))]
%
% INPUTS:
%   Z  -   mz-by-(k+3) real: i/o data
%   p  -   mp-by-(k+2) integer>=0: determines which monomials to use
% 
% OUTPUTS:
%   f   -  scalar msspoly in z=[v;w]: v=[v0;v1], w=[w0,...,wk]
%   h   -  scalar msspoly in v1
%   r   -  scalar msspoly in z
%   U   -  column of msspoly monomials in variables z,d=[d0;d1]
%   Q   -  symmetric positive semidefinite matrix
% 
% r,h,U,Q certify robustness of (*) via the polynomial identity
%   r+d(1)*(f+diff(f,v,d))+h(v(2))*d(2)^2-h(v(1))*d(1)^2-d(1)^2=U'*Q*U.
% Robust fidelity measure sum_i r(Z(i,:)) is minimized.
%
% p  defines U=[U0;U1*d(1);U2*d(2)], where U0,U1,U2
%    are vectors of monomials in z:
%
% U0 has monomials (v0^a0)*(v1^a1)*(w0^b0)*...*(wk^bk), where
%     either a0<=p(1), a1=0, and b0/p(3)+b1/p(4)+...+bk/p(k+3)<=1     
%     or 2*a0<=p(1),2*a1<=p(1), and b0/p(3)+b1/p(4)+...+bk/p(k+3)<=1/2
% U1 has monomials (v0^a0)*(v1^a1)*(w0^b0)*...*(wk^bk), where
%     either a0<p(1), a1=0, and b0/p(3)+b1/p(4)+...+bk/p(k+3)<=1     
%     or 2*a0<p(1),2*a1<=p(1), and b0/p(3)+b1/p(4)+...+bk/p(k+3)<=1/2
% U2 has monomials v1^a1, where a1<p(1)

if nargin<2, error('2 inputs required'); end
if ~isa(Z,'double'), error('input 1 not a "double"'); end
[mz,nz]=size(Z);
if mz<1, error('input 1 is empty'); end
if ~mint_isint(p), error('input 2 not integer'); end
[mp,np]=size(p);
if mp~=1, error('input 2 is not a row vector'); end
if np~=nz-1, error('inputs 1,2 not compatible'); end
k=np-3;
if k<0, error('input 2 is too short'); end

v=msspoly('v',2);          % absract variables
w=msspoly('w',k+1); 
z=[v;w];

Kw=[zeros(1,k+1);diag(p(2:k+2))];
K0=[repmat(0,k+2,1) Kw;repmat(p(1),k+2,1) Kw]
p0=mint_ch();
U0=recomp(z,p0);
m0=size(U0,1);
p1=mint_ch([repmat(0,k+2,1) Kw;repmat(p(1)-1,k+2,1) Kw]);
U1=recomp(z,p1);
m1=size(U1,1);

fprintf('\nnlid_miso0: pre-computing.');tic;
C0=zeros(m0); 
for i=1:mz, 
    u0=double(subs(U0,z,Z(i,:)'));
    C0=C0+u0*u0'; 
    if mod(i,100)==0, fprintf('.'); end
end
toc
pr=mssprog;
[pr,Q]=new(pr,m0+m1,'psd','Q');
Q11=Q(1:m0,1:m0);
Q12=Q(1:m0,m0+1:m0+m1);
Q22=Q(m0+1:m0+m1,m0+1:m0+m1);
r=U0'*Q11*U0;
f=U0'*Q12*U1;
e=U1'*Q22*U1+1-2*diff(f,v);
pr=eq(pr,e);
pr=sedumi(pr,trace(Q11*C0));
f=get(pr,f);
r=get(pr,r);
U=[U0;U1];
Q=double(pr,Q);
fprintf('nlid_miso: upper bound %e, certified: %f>0\n', ...
    trace(double(pr,Q11)*C0),min(eig(Q)))
