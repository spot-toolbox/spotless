function y=sim_findzero(p,t,k)
% function y=sim_findzero(p,t,k)
%
% solve equation f(y)=0 where f:R->R has positive derivative
%
% INPUTS:
%   p  -  function handle: p(y)=[f(y) f'(y)], assumed f'>0
%   t  -  1-by-1 real:     initial guess (default t=0)
%   k  -  1-by-1 integer:  number of iterations (default k=5)
% OUTPUT:
%   y  -  1-by-1 real:     a solution of  f(y)=0

if isa(p,'double'),           % to readily accept polynomials 
    p=p(:)';
    p=@(x)[polyval(p,x) polyval(p(1:length(p)-1).*(length(p)-1:-1:1),x)];
end
if nargin<2, t=0; end
if nargin<3, k=5; end

y=t;  
zy=p(y);
for i=1:k,
    u=y-zy(1)/(zy(2)+(zy(2)==0));
    zu=p(u);
    if zy(1)*zu(1)<=0, break; end
    y=u;
    zy=zu;
end
if (zu(1)==0)||(i==k), y=u; return; end
if y<u,
    a=y; b=u;
else
    a=u; b=y;
end
y=0.5*(a+b);   % a<y<b ALWAYS
for j=i:k,
    z=p(y);           % calculate value/derivative
    if z(1)==0,       % update bounds
        return
    elseif z(1)>0,
        b=y;
    else
        a=y;
    end
    y=y-z(1)/(z(2)+(z(2)==0));    % update y
    if (y<=a)||(y>=b), y=0.5*(a+b); end
end