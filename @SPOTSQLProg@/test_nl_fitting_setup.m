%
%
% Compare primal and dual standard form solutions for a fitting problem:
%
%  [ r_t  0        e_t'
%    0    E-C'C  F_t'
%    e_t  F_T      E   ]
%
%  e_t = Ex(t+1) - f(x(t),u(t))
%  F_t = df/dx(x(t),u(t))
%

m = 1;
n = 6;

A = randn(n,n);
A = A*0.9/max(abs(eig(A)));

B = randn(n,1)/5;

C = randn(1,n)/5;
Q = C'*C;

R = dlyap(A',Q);

N = 1e5;
us = 0.5*randn(1,N);
xs = zeros(n,N+1);
for i = 1:N
    xs(:,i+1) = A*xs(:,i)  + B*us(:,i);
end
xs = xs + 0.001*randn(size(xs));

U = us(:,1:N).';
X = xs(:,1:N).';
V = xs(:,2:N+1).';


x = msspoly('x',n);
v = msspoly('v',n);
u = msspoly('u',m);
