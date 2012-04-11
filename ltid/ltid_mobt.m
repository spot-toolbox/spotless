function H=ltid_mobt(G,r)
% function H=ltid_mobt(G,r)
%
% H(z)=G((z-r)/(1-rz))

h=sqrt(1-r^2);
[A,B,C,D]=ssdata(G);
I=eye(size(A,1));
X=inv(I+r*A);
c=C*X;
d=D-r*(c*B);
c=c*h;
b=X*B*h;
a=(r*I+A)*X;
H=ss(a,b,c,d,-1);
