function [h,r]=dsgn_dtsf(d,u,w)
% function [h,r]=dsgn_dtsf(d,u,w)
%
% optimized low-pass filter H=H(z) with frequency response 
%    H(exp(jt))=cos(t*d)*h, 
% minimizing r=max(abs(cos(w*d)*h)) subject to min(cos(u*d)*h)>1
%
% INPUTS:
%    d  -  1-by-k integer
%    u  -  n-by-1 real
%    w  -  N-by-1 real

if nargin<3, error('three inputs required'); end
k=size(d,2);
if size(d,1)~=1, error('input 1 is not a row vector'); end
u=real(u(:));
w=real(w(:));
n=length(u);
N=length(w);

h=msspoly('h',[k 1]);           % coefficients of h
r=msspoly('r');                 % minimization objective
x=msspoly('x',n+2*N);           % slack variables 
hu=cos(u*d)*h;                  % samples at u
hw=cos(w*d)*h;                  % samples at w

pr=mssprog;                     % initialize optimization setup
pr.free=h;                      % register h as free
pr.free=r;                      % register r as free
pr.pos=x;                       % register x as positive (point-wise)
pr.eq=hu-1-x(1:n);              % enforce hu>1
pr.eq=r-hw-x(n+1:n+N);          % enforce r>hw
pr.eq=r+hw-x(n+N+1:n+2*N);      % enforce r>-hw

pr.sedumi=r;                    % minimize r

h=pr({h});
r=pr({r});