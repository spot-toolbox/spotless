function [x,i]=ltid_sfls(a,b,h,r,k)
% function [x,i]=ltid_sfls(a,b,h,r,k)
%
% functional least squares: tries to minimize|Ax-B| subject to |x|<r|B|
%
% INPUTS:
%   a  -  function handle a=@(x)(A'*(A*x))
%   b  -  b=A'*B
%   h  -  |B|
%   r  -  positive real scalar (default r=1e6)
%   k  -  positive integer (maximal number of iterations, default 200)
%
% OUTPUTS:
%   x  -  n-by-1  (|x|<r|B| and 4r|A'(Ax-B)|<|B|  when converged)
%   i  -  integer scalar (numit for |x|<r|B|, -numit otherwise)

if nargin<3, error('3 inputs required'); end
if nargin<4, r=1e6; end
if nargin<5, k=200; end
if ~isa(a,'function_handle'), error('input 1 not a function handle'); end
if ~isa(b,'double'), error('input 2 not a double'); end
[n,nb]=size(b);
if nb~=1, error('input 3 not a column'); end

x=zeros(n,1);               % current guess
i=0;                        % current step
if h<=0, return; end        
b=b/h;                      % normalize
t=1/(4*r);                  % absolute tolerance after normalization
u=zeros(n,k);               % orthogonalized basis in R(A), columnwise
v=zeros(n,k);               % a*v=u

w=a(b);
rw=norm(w);
w=w/rw;
q=b/rw;
for i=1:k,
    u(:,i)=w;
    v(:,i)=q;
    z=(b'*u(:,1:i))';
    b=b-u(:,1:i)*z;
    x=x+v(:,1:i)*z;
    rb=norm(b);
    if rb<t, break; end
    rx=norm(x);
    if rx>r, i=-i; break; end
    w=a(b);
    z=(w'*u(:,1:i))';
    w=w-u(:,1:i)*z;
    q=b-v(:,1:i)*z;
    rw=norm(w);
    if rw<=0, break; end
    w=w/rw;
    q=q/rw;   
end
x=h*x;
    
