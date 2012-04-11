function [g,y]=mint_gcd(x)
% function [g,y]=mint_gcd(x)
%
% For a non-zero double integer m-by-1 vector x, g=y*x is the 
% minimal (over all non-integer z) non-zero value of abs(z*x),
% achieved, in particular, at z=y. For x=0, g=0 and  y=0.

if nargin<1, error('one input required'); end
if ~mint_isint(x), error('input is not a double integer'); end
[m,n]=size(x);  
if n~=1, error('input not a column vector'); end
if m==1, 
    g=abs(x);
    y=sign(x+0.5); 
elseif m==2,
    [g,c,d]=gcd(x(1),x(2));
    y=[c d];
else
    if all(x==0),
        g=0;
        y=x';
    else
        [g0,y0]=mint_gcd(x(2:m));
        [g,c,d]=gcd(g0,x(1));
        y=[d c*y0];
        [r,i]=max(1./(abs(x)-0.5));
        y=mod(y,abs(x(i)));
        y(i)=y(i)+(g-y*x)/x(i);
    end
end
