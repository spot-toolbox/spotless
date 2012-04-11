function [f,r,U0,U1,Q]=nlid_miso0(Z,p)
% function [f,r,U0,U1,Q]=nlid_miso0(Z,p)
%
% Non-Linear System Identification for 
% Multiple Input Single Output systems (zero order version)
%
% Generate implicit polynomial MISO model (input w, output v)
%       f(v(t),w0(t),...,wk(t))=0           (*)
% from data Z(i,:)=[v(t(i)) w0(t(i)) ... wk(t(i))]
%
% INPUTS:
%   Z  -   mz-by-(k+2) real: i/o data
%   p  -   mp-by-(k+2) integer>=0: determines which monomials to use
% 
% OUTPUTS:
%   f   -  scalar msspoly in z=[v;w]: v=msspoly('v'), w=msspoly('w',k+1)
%   r   -  scalar msspoly in z
%   U0  -  column of msspoly monomials in z
%   U1  -  column of msspoly monomials in z
%   Q   -  symmetric positive semidefinite matrix
% 
% r,U0,U1,Q certify robustness of (*) via the polynomial identity
%   r+d*(f+diff(f,v,d))-d^2=U'*Q*U,  where U=[U0;U1*d], d=msspoly('d')
% Robust fidelity measure sum_i r(Z(i,:)) is minimized.
%
% p  defines U0,U1 as vectors of monomials in z:
%
% U0 has monomials from  (v^a0)*(w0^b0)*...*(wk^bk),
%     where a0<=p(1) and b0/p(2)+b1/p(3)+...+bk/p(k+2)<=1
% U1 has monomials of the form  (v^a0)*(w0^b0)*...*(wk^bk),
%     where a0<p(1) and b0/p(2)+b1/p(3)+...+bk/p(k+2)<=1

if nargin<2, error('2 inputs required'); end
if ~isa(Z,'double'), error('input 1 not a "double"'); end
[mz,nz]=size(Z);
if mz<1, error('input 1 is empty'); end
if ~mint_isint(p), error('input 2 not integer'); end
if any(p(:)<1), error('input 2 not positive'); end
[mp,np]=size(p);
if mp~=1, error('input 2 is not a row vector'); end
if np~=nz, error('inputs 1,2 not compatible'); end
k=np-2;
if k<0, error('input 2 istoo short'); end

v=msspoly('v');          % absract variables
w=msspoly('w',k+1); 
z=[v;w];

Kw=[zeros(1,k+1);diag(p(2:k+2))];
p0=mint_ch([repmat(0,k+2,1) Kw;repmat(p(1),k+2,1) Kw]);
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
Q=double(pr,Q);
fprintf('nlid_miso: upper bound %e, certified: %f>0\n', ...
    trace(double(pr,Q11)*C0),min(eig(Q)))

