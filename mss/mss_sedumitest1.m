function mss_sedumitest1
% function mss_sedumitest1
%
% minimize x1 subject to [x1 x2;x3 x4]>0, x1+x2+x3+x4=0

A=[1 1 1 1];
B=0;
C=[1 0 0 0]';
K.s=2;
[x,y,info]=sedumi(A,B,C,K)