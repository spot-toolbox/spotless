function [a,b,P,Q]=nlid_miso_old(vv,ww,F,R,U,Z)
% function [a,b,P,Q]=nlid_miso_old(vv,ww,F,R,U,Z)
%
% Non-Linear System Identification for Multiple Input Single Output systems
%
% INPUTS:
%   vv  -  mt-by-r real (rows are output stack samples)
%   ww  -  mt-by-k real (rows are input samples)
%   F   -  mf-by-1 polynomial in v=msspoly('v',r), w=msspoly('w',k)
%   R   -  mr-by-1 polynomial in v,w
%   U   -  mu-by-1 polynomial in v,w,d=msspoly('d',r), affine in d
%   Z   -  mz-by-1 polynomial in v(2:r),d(2:r), linear in d
% 
% OUTPUTS:
%   a   -  1-by-mf real
%   b   -  1-by-mr real
%   P   -  mz-by-mz real symmetric positive semidefinite
%   Q   -  mu-by-mu real symmetric positive semidefinite
% 
% Define polynomials f=a*F, r=b*R, h0=Z'*P*Z, g=U'*Q*U
% such that the model MISO model defined by equations
%   f(v(t),v(t-1),...,v(t-r+1),w(t))=0
% has its robust local incremental stability certified by the identity
%   d(1)*(f+diff(f,v,d))-d(1)^2+r(v,w)+h0-h1=g, 
% where h1 is the result of substituting 
% [v(1:r-1);d(1:r-1)] for [v(2:r);d(2:r)] in h0.
% The robust fidelity measure sum_i r(vv(i,:),w(i,:)) is minimized


if nargin<5, error('5 inputs required'); end
if ~isa(vv,'double'), error('1st input not a "double"'); end
if isempty(vv), error('1st input is empty'); end
if ~isa(ww,'double'), error('2nd input not a "double"'); end
if isempty(ww), error('2nd input is empty'); end
if ~all(isfinite(vv(:))), error('1st input not finite'); end
if ~all(isfinite(ww(:))), error('2nd input not finite'); end
[T,r]=size(vv); 
[mww,k]=size(ww); 
if T~=mww, error('1st and 2nd inputs are incompatible'); end
if ~isa(F,'msspoly'), error('3rd argument not a mss polynomial'); end
if ~isa(R,'msspoly'), error('4th argument not a mss polynomial'); end
if ~isa(U,'msspoly'), error('5th argument not a mss polynomial'); end
[mf,nf]=size(F); 
[mr,nr]=size(R); 
[mu,nu]=size(U); 
if nf~=1, error('3rd input must be a column'); end
if nr~=1, error('4th input must be a column'); end
if nu~=1, error('5th input must be a column'); end
if r>1,
    if nargin<6, error('6 inputs required for models with memory'); end
    if ~isa(Z,'msspoly'), error('6th argument not an mss polynomial'); end
    [mz,nz]=size(Z); 
    if nz~=1, error('6th input must be a column'); end
end    
v=msspoly('v',r);  
w=msspoly('w',k); 
d=msspoly('d',r);
vw=[v;w];
vwd=[v;w;d];
if ~isfunction(F,vw), error('illegal variables in the 3rd input'); end
if ~isfunction(R,vw), error('illegal variables in the 4th input'); end
if ~isfunction(U,vwd), error('illegal variables in the 5th input'); end
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

