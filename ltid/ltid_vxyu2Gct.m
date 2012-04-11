function G=ltid_vxyu2Gct(V,X,Y,U,h,o)
% function G=ltid_vxyu2Gct(V,X,Y,U,h,o)
%
% convert input/state/output samples to stable CT SS LTI model
%
% INPUTS:
%  V  -  n-by-N real: samples of dx/dt (
%  X  -  n-by-N real: samples of x
%  Y  -  k-by-N real: samples of y
%  U  -  m-by-N real: samples of u
%  h  -  positive scalar (tolerance for strict positivity, default 1e-4)
%  o  -  display flag (default o=0, no display)
%
% OUTPUTS:
%  G  -  stable CT LTI state space model
%
% fits equations E*dx/dt+Px+Bu=0, C*dx/dt+Du=y to data
% by minimizing tr(Q[V;X;Y;U][V;X;Y;U]') subject to P=P'>0 and
% positive definiteness of 
%   f(v,x,y,u)+2w'(Ew+Ev+Px+Bu)+d'd+2d'(Cw+Cx+Du-y)

if nargin<4, error('4 inputs required'); end   
if nargin<5, h=1e-4; end
if nargin<6, o=0; end
[n,N]=size(V);
if ~isequal([n,N],size(X)), error('inputs 1,2 incompatible'); end
[k,nn]=size(Y);
if nn~=N, error('inputs 1,3 incompatible'); end
[m,nn]=size(U);
if nn~=N, error('inputs 1,4 incompatible'); end
n2=2*n;                            % dimension of [v;x]
n3=n2+k;                           % dimension of [v;x;y]
n4=n3+m;                           % dimension of [v;x;y;u]
n5=n4+n;                           % dimension of [v;x;y;u;w]
n6=n5+k;                           % dimension of [v;x;y;u;w;d]
W=[V;X;Y;U]; 
W=W*W';                            % correlation matrix
pr=mssprog;                        % initialize mssprog 
Q=msspoly('Q',nchoosek(n6+1,2));   % n6-by-n6 Q=Q'>0
pr.psd=Q;
Q=mss_v2s(Q);  
P=msspoly('P',nchoosek(n+1,2));    % n-by-n P=P'>0
pr.psd=P;
P=mss_v2s(P);
% Q(n4+1:n6,:)=[ E      P      zeros(n,k)  B   E+E'-h*eye(n)    C'  ; ...
%              [ C  zeros(k,n)   -eye(k)   D        C         eye(k)]
E=Q(n4+1:n5,1:n);
pr.eq=Q(n4+1:n5,n+1:n2)-P;
pr.eq=Q(n4+1:n5,n2+1:n3);
B=Q(n4+1:n5,n3+1:n4);
pr.eq=mss_s2v(E+E'-Q(n4+1:n5,n4+1:n5)-h*eye(n));
C=Q(n5+1:n6,1:n);
pr.eq=Q(n5+1:n6,n+1:n2);
pr.eq=Q(n5+1:n6,n2+1:n3)+eye(k);
D=Q(n5+1:n6,n3+1:n4);
pr.eq=Q(n5+1:n6,n4+1:n5)-C;
pr.eq=mss_s2v(Q(n5+1:n6,n5+1:n6)-eye(k));
y=trace(W*Q(1:n4,1:n4));
pr=sedumi(pr,y,o);
Ei=inv(pr({E}));
P=pr({P});
B=pr({B});
C=pr({C});
D=pr({D});
A=-Ei*P;
B=-Ei*B;
G=ss(A,B,C*A,D+C*B);