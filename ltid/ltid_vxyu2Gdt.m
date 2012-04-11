function G=ltid_vxyu2Gdt(V,X,Y,U,h,o)
% function G=ltid_vxyu2Gdt(V,X,Y,U,h,o)
%
% convert input/state/output samples to stable DT SS LTI model
%
% INPUTS:
%  V  -  n-by-N real: samples of x[t+1] 
%  X  -  n-by-N real: samples of x[t]
%  Y  -  k-by-N real: samples of y[t]
%  U  -  m-by-N real: samples of u[t]
%  h  -  positive scalar (tolerance for strict positivity, default 1e-4)
%  o  -  display flag (default o=0, no display)
%
% OUTPUTS:
%  G  -  stable DT LTI state space model
%
% fits equations Ex[t+1]+Ax[t]+Bu[t]=0, Cx[t]+Du[t]=y[t] to data
% by minimizing tr(Q[V;X;Y;U][V;X;Y;U]') subject to H=H'>0 and
% positive definiteness of 
%   f(v,x,y,u)+2w'(Ew+Ae+Ev+Ax+Bu)+e'He-w'Hw+d'd+2d'(Ce+Cx+Du-y)

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
n7=n6+n;                           % dimension of [v;x;y;u;w;d;e]
W=[V;X;Y;U]; 
W=W*W';                            % correlation matrix
pr=mssprog;                        % initialize mssprog 
Q=msspoly('Q',nchoosek(n7+1,2));   % n7-by-n7 Q=Q'>0
pr.psd=Q;
Q=mss_v2s(Q);  
% Q(n4+1:n7,:)=[   E      A    O(n,k)    B    E+E'-H-h*I(n) O(n,k) A ; ...
%              [ O(k,n)   C    -I(k)     D       O(k,n)      I(k)  C ; ...
%              [  O(n)   O(n)  O(n,k)  O(n,m)      A'         C'   H ]
E=Q(n4+1:n5,1:n);                                     % (5,1)
A=Q(n4+1:n5,n+1:n2);                                  % (5,2)
pr.eq=Q(n4+1:n5,n2+1:n3);                             % (5,3)
B=Q(n4+1:n5,n3+1:n4);                                 % (5,4)
pr.eq=Q(n5+1:n6,1:n);                                 % (6,1)
C=Q(n5+1:n6,n+1:n2);                                  % (6,2)
pr.eq=Q(n5+1:n6,n2+1:n3)+eye(k);                      % (6,3)
D=Q(n5+1:n6,n3+1:n4);                                 % (6,4)
pr.eq=Q(n5+1:n6,n4+1:n5);                             % (6,5)
pr.eq=mss_s2v(Q(n5+1:n6,n5+1:n6)-eye(k));             % (6,6)
H=Q(n6+1:n7,n6+1:n7);                                 % (7,7)
pr.eq=mss_s2v(E+E'-H-Q(n4+1:n5,n4+1:n5)-h*eye(n));    % (5,5)
pr.eq=Q(n6+1:n7,1:n4);                                % (7,1:4)
pr.eq=Q(n6+1:n7,n4+1:n5)-A';                          % (7,5)
pr.eq=Q(n6+1:n7,n5+1:n6)-C';                          % (7,6)
y=trace(W*Q(1:n4,1:n4));
pr=sedumi(pr,y,o);
Ei=inv(pr({E}));
A=pr({A});
B=pr({B});
C=pr({C});
D=pr({D});
G=ss(-Ei*A,-Ei*B,C,D,-1);