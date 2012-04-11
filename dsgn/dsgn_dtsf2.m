function [h,r]=dsgn_dtsf2(T,m,N)
% function [h,r]=dsgn_dtsf2(T,m,N)
%
% optimize coefficients h=[hm;0;...;h2;0;h1;1;h1;0;h2;0;h3;0;...;hm] 
% of a "smoothing filter": 
% i.e. minimize max{|H(t)|: pi-T<t<pi} for H(t)=cos((-2*m+1:2*m-1)*t)*h
%
% for a DT signal x[] with zero odd samples, y=conv(x,h) has
% same even samples, and minimal spectrum in the [pi-T,pi] range

if nargin<1, T=pi/4; end
if (T<=0)||(T>=pi/2), error('T out of bounds'); end
if nargin<2, m=4; end
if (m~=round(m))||(m<1), error('m is not a positive integer'); end
if nargin<3, N=300; end
cs=cos(linspace(pi-T,pi,N)'*(1:2:(2*m-1)));

h=msspoly('h',[m 1]);           % coefficients of h
r=msspoly('r');                 % minimization objective
e=cs*h+0.5;                     % samples of H+0.5
x=msspoly('x',2*N);             % slack variables 

pr=mssprog;                     % initialize optimization setup
pr.free=h;                      % register h as free
pr.free=r;                      % register r as free
pr.pos=x;                       % register x as positive (point-wise)
pr.eq=e+r-x(1:N);               % enforce H+0.5>-r
pr.eq=r-e-x(N+1:2*N);           % enforce H+0.5<r

pr.sedumi=r;                    % minimize r

h0=pr({h});
h=zeros(4*m-1,1);
h(1:2:2*m-1)=h0(m:-1:1);
h(2*m)=1;
h(2*m+1:2:4*m-1)=h0;
r=2*pr({r});