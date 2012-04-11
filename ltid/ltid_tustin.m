function Gc=ltid_tustin(G,f0)
% function Gc=ltid_tustin(G,f0)
%
% Tustin transform with z=(1+s/f0)/(1-s/f0)

if nargin<2, error('2 inputs required'); end
[a,b,c,d]=ssdata(G);
I=eye(size(a,1));
e=inv(I+a);
D=d-c*e*b;
A=f0*(a-I)*e;
B=(2*f0)*e*b;
C=c*e;
Gc=ss(A,B,C,D);