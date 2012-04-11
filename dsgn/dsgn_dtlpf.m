function h=dsgn_dtlpf(m,N,T)
% function h=dsgn_dtlpf(m,N,T)
%
% optimize  H(t)=sin((-m:m)*t)*h (h=[-hm;...;-h2;-h1;0;h1;h2;...;hm])
% with two arguments:
%    minimize sum of squares of H(t)-t over t=linspace(0,pi,N) 
% with three arguments (0<T<pi/2):
%    minimize max of |H(t)-t| over t=linspace(0,pi-T,N) 

if nargin<1, m=200; end
if nargin<2, N=3000; end

if nargin<3, 
    t=linspace(0,pi,N)';
    sn=sin(t*(1:m));
    h0=0.5*((sn'*sn)\(sn'*t));
else
    t=linspace(0,pi-T,N)';
    sn=sin(t*(1:m));
    h=msspoly('h',[m 1]);           % coefficients of h
    r=msspoly('r');                 % minimization objective
    e=sn*h-0.5*t;                   % samples of (H(t)-t)/2
    x=msspoly('x',2*N);             % slack variables 
    pr=mssprog;                     % initialize optimization setup
    pr.free=h;                      % register h as free
    pr.free=r;                      % register r as free
    pr.pos=x;                       % register x as positive (point-wise)
    pr.eq=e+r-x(1:N);               % enforce (H-t)/2>-r
    pr.eq=r-e-x(N+1:2*N);           % enforce (H-t)/2<r
    pr.sedumi=r;                    % minimize r
    h0=pr({h});
end
h=zeros(2*m+1,1);
h(1:m)=-h0(m:-1:1);
h(m+1)=0;
h(m+2:2*m+1)=h0;