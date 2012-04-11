function y=nlid_iosim(f,u,y0,r)
% function y=nlid_iosim(f,u,y0,r)
%
% calculate response y=[y(1);...;y(mu)]  (y(t)'s are scalars)
% to input u=[w{1};...;w{mu}] (w{t}'s are 1-by-nu)
% of the MISO system f(y(t),y(t-1),...,y(t-my),w(t))=0, 
% with initial conditions y0=[y(1);...;y(my)],
% using r (default r=1) iterations of the Newton method 
% to solve f(y(t),y(t-1),...,y(t-my),w{t})=0 for y(t) at each step
%
% f must be mss polynomial of v=msspoly('v',my+1) and w=msspoly('w',nu)
% it is expected that the polynomial diff(f,v(1)) has no roots at all

if nargin<3, error('2 inputs required'); end
if ~isa(u,'double'), error('2nd input not a "double"'); end
if isempty(u), error('2nd input is empty'); end
if ~isa(y0,'double'), error('3rd input not a "double"'); end
if isempty(y0), error('3rd input is empty'); end
[mu,nu]=size(u);
[my,ny]=size(y0);
if ny~=1, error('3rd input is not a column'); end
if nargin<4, r=1; end; 
if ~isa(r,'double'), error('4th input not a "double"'); end
if ~isscalar(r), error('4th input not a scalar'); end
if r~=max(1,round(r)), error('4th input not integer>0'); end
v=msspoly('v',my+1);
w=msspoly('w',nu);
vw=[v;w];
vw0=[v(2:my+1);w];
if ~isfunction(f,vw), error('illegal variables in the 1st input'); end
% q=diff(f,v(1));
y=zeros(mu+my,1);
y(1:my)=y0;
fprintf('nlid_iosim.')
for t=1:mu,
    xt=[y(t+my-1:-1:t);u(t,:)'];
    F=subs(f,vw0,xt);
    y(t+my)=newton(F,v(1),r,y(t+my-1));
    if mod(t,100)==0, fprintf('.'); end
end
fprintf('done.\n')

    