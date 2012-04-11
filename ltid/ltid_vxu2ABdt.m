function [A,B]=ltid_vxu2ABdt(V,X,U,h,o)
% function [A,B]=ltid_vxu2ABdt(V,X,U,h,o)
%
% convert input/state/output samples to stable x(t+1)=Ax(t)+Bu(t)
%
% INPUTS:
%  V  -  n-by-N real: samples of x[t+1] 
%  X  -  n-by-N real: samples of x[t]
%  U  -  m-by-N real: samples of u[t]
%  h  -  positive scalar (tolerance for strict positivity, default 1e-4)
%  o  -  display flag (default o=0, no display)
%
% OUTPUTS:
%  A  -  n-by-n real Schur matrix
%  B  -  n-by-m real
%
% fits equations L*[x(t+1);x(t);u(t)]=Ex(t+1)+ax(t)+bu(t)=0 (a=-EA,b=-EB)
% to data by minimizing tr(L[V;X;U][V;X;U]'L') subject to 
%     H=H'>0,   [H a';a E+E'-H]>0,    trace(E)=1

if nargin<3, error('3 inputs required'); end   
if nargin<4, h=1e-4; end
if nargin<5, o=0; end
[n,N]=size(V);
if ~isequal([n,N],size(X)), error('inputs 1,2 incompatible'); end
[m,nn]=size(U);
if nn~=N, error('inputs 1,3 incompatible'); end

%[U,S,V]=svd([V;X;U],'econ'); 
%W=U*S*U';                         % square root of [V;X;U]*[V' X' U']
W=[V;X;U];
clear V X U
[Q,R]=qr(W',0);
W=R';
clear Q R

pr=mssprog;                       % initialize mssprog 
Q=msspoly('Q',nchoosek(2*n+1,2));   % 2n-by-2n Q=[H A';A E+E'-H]
pr.psd=Q;
Q=mss_v2s(Q); 
H=Q(1:n,1:n);
A=Q(n+1:2*n,1:n);
E=msspoly('E',n^2);
pr.free=E;
E=reshape(E,n,n);
pr.eq=1-trace(E);
B=msspoly('B',n*m);
pr.free=B;
B=reshape(B,n,m);
pr.eq=h*eye(n)+Q(n+1:2*n,n+1:2*n)+H-E-E';
nx=n*(2*n+m)+1;
x=msspoly('x',nx);
pr.lor=x;
pr.eq=x(2:nx)-reshape([E A B]*W,nx-1,1);

pr=sedumi(pr,x(1),o);
Ei=inv(pr({E}));
A=-Ei*pr({A});
B=-Ei*pr({B});