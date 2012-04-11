function x=ltid_Linf(A,B)
% function x=ltid_Linf(A,B)
%
% given  m-by-n real matrix A with m>n, and a real m-by-1 vector B,
% x minimizes the L-infinity norm of Ax-B

[m,n]=size(A);                % dimensions of A
x=msspoly('x',[n 1]);         % the actual decision vector
y=msspoly('y');               % objective variable
z=msspoly('z',[2*m 1]);       % slack variables
e=A*x-B;                      % fitting error 
pr=mssprog;                   % initialize mss program
pr.free=x;                    % x is free
pr.free=y;                    % y is free
pr.pos=z;                     % z>0
pr.eq=y+e-z(1:m);             % y>-e 
pr.eq=y-e-z(m+1:2*m);         % y>e
pr.sedumi=y;                  % minimize y
x=pr({x});                    % get the answer as a "double"