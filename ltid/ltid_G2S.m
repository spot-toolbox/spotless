function S=ltid_G2S(G)
% function S=ltid_G2S(G)
%
% D=(I+G)^{-1}(I-G)

[A,B,C,D]=ssdata(G);
I=eye(size(D,1));
d=inv(I+D);
c=-d*C;
d=2*d;
b=B*d;
a=A+B*c;
S=minreal(ss(a,b,c,d-I));