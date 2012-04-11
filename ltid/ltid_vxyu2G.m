function [G,obj,E,F,L,C,D,P]=ltid_vxyu2G(V,X,Y,U,T,h,o)
% function [G,obj,E,F,L,C,D,P]=ltid_vxyu2G(V,X,Y,U,T,h,o)
%
% convert input/state/output samples to stable SS LTI model
%
% INPUTS:
%  V  -  n-by-N real: V(:,i)=(x[ti+T]-x[ti])/T  (T=0 admissible)
%  X  -  n-by-N real: X(:,i)=x[ti]
%  Y  -  k-by-N real: Y(:,i)=y[ti]
%  U  -  m-by-N real: U(:,i)=u[ti]
%  T  -  positive real: sampling time (T=0 means a CT data), default T=1
%  h  -  positive scalar (tolerance for strict positivity, default h=1e-4)
%  o  -  algorithm/display flag (default o=0, no display)
%
% OUTPUTS:
%  G  -  stable LTI state space model (CT when T=0, DT otherwise)
%
% fits LTI state space equations 
%   E(x[t+T]-x[t])=T*(Fx[t]+Lu[t]), y[t]=Cx[t]+Du[t]    (when T>0)
%           Ex'[t]=Fx[t]+Lu[t],     y[t]=Cx[t]+Du[t]    (when T=0)
% to (V,X,Y,U) data by minimizing obj=tr(R*([V;X;Y;U][V;X;Y;U]')) subj. to
% r+e'Pe+|q|^2+2z'(Cz+a)+z'Pz+2z'(Ee-Pe+.5h(Fe+b))+w'Pw+2w'(Fe+Pe+b) >> 0
% for all v,x,y,u,e,z,w,q, for a=Fx+Lu-Ev,b=Cx+Du-y,r=[v;x;y;u]'R[v;x;y;u]


if nargin<4, error('4 inputs required'); end 
if nargin<5, T=1; end
if nargin<6, h=1e-4; end
if nargin<7, o=1; end
[n,N]=size(V);
if ~isequal([n,N],size(X)), error('inputs 1,2 incompatible'); end
[k,nn]=size(Y);
if nn~=N, error('inputs 1,3 incompatible'); end
[m,nn]=size(U);
if nn~=N, error('inputs 1,4 incompatible'); end
T=real(T(1));
if T<0, error('input 5 is negative'); end
n2=2*n;                            % dimension of [v;x]
n3=n2+m;                           % dimension of [v;x;u]
n4=n3+k;                           % dimension of [v;x;u;y]
n5=n4+n;                           % dimension of [v;x;u;y;e]
n6=n5+n;                           % dimension of [v;x;u;y;e;z]
n7=n6+n;                           % dimension of [v;x;u;y;e;z;w]
n8=n7+k;
W=[V;X;U;Y]; 
W=W*W';                            % correlation matrix
wmx=max(abs(W(:)));
W=W/wmx;
pr=mssprog;                        % initialize mssprog 
Q=msspoly('Q',nchoosek(n8+1,2));   % n8-by-n8 Q=Q'>0
pr.psd=Q;
Q=mss_v2s(Q);  
%  Q=[R H';H S],
%  H=[  O(n)    O(n)  O(n,m)    O(n,k) ] 
%    [ -ET/2   FT/2    LT/2     O(n,k) ] 
%    [   -E      F       L      O(n,k) ]   
%    [O(k,n)   C      D      -I(k) ]    
%  S=[    P      E'-P+.5F'T   F'+P     C'   ]
%    [ E-P+.5FT      P        O(n)   O(n,k) ]
%    [   F+P        O(n)       P     O(n,k) ]
%    [    C        O(k,n)    O(k,n)   I(k)  ]
E=-Q(n6+1:n7,1:n);                        % matrices of interest
F=Q(n6+1:n7,n+1:n2);
FT=(T/2)*F;
L=Q(n6+1:n7,n2+1:n3);
P=Q(n6+1:n7,n6+1:n7);
Pv=mss_s2v(P);
C=Q(n7+1:n8,n+1:n2);
D=Q(n7+1:n8,n2+1:n3);
pr.eq=Q(n4+1:n5,1:n4);                    % zero entries
pr.eq=Q(n7+1:n8,1:n);
pr.eq=Q(n5+1:n7,n3+1:n4);
pr.eq=Q(n6+1:n8,n5+1:n6);
pr.eq=Q(n7+1:n8,n6+1:n7);
pr.eq=mss_s2v(Q(n4+1:n5,n4+1:n5))-Pv;     % other P
pr.eq=mss_s2v(Q(n5+1:n6,n5+1:n6))-Pv;
pr.eq=Q(n7+1:n8,n3+1:n4)+eye(k);          % identity matrices
pr.eq=Q(n7+1:n8,n7+1:n8)-eye(k);
pr.eq=(T/2)*Q(n6+1:n7,1:n3)-Q(n5+1:n6,1:n3);
pr.eq=Q(n7+1:n8,n4+1:n5)-C;
pr.eq=E-P+FT-Q(n5+1:n6,n4+1:n5);
pr.eq=F+P-Q(n6+1:n7,n4+1:n5);
obj=trace(W*Q(1:n4,1:n4));
pr=sedumi(pr,obj,o>0);
E=pr({E});
F=pr({F});
L=pr({L});
C=pr({C});
D=pr({D});
P=pr({P});
Ei=inv(E);
obj=wmx*pr({obj});
if T==0,
    G=ss(Ei*F,Ei*L,C,D);
else
    Ei=T*Ei;
    G=ss(eye(n)+Ei*F,Ei*L,C,D,-1);
end