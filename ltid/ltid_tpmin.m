function [y,t]=ltid_tpmin(a)
% function [y,t]=ltid_tpmin(a)
%
% INPUT:
%   a  -  m-by-1 real
%
% OUTPUT:
%   y  -  minimum of cos((0:m-1)*t)*a over real t
%   t  -  argument of minimum of cos((0:m-1)*t)*a over real t

if nargin<1, error('1 input required'); end
if ~isa(a,'double'), error('input not a double'); end
if ~isreal(a), error('input not real'); end
[m,n]=size(a); 
if (n~=1)||(m<1), error('input 1 not a column vector'); end
if m==1,                   % the constant polynomial case
    y=a(1);
    t=0;
    return
end
ad=(1:m-1)'.*a(2:m);
w=angle(roots([-ad(m-1:-1:1);0;ad]));
[y,i]=min(cos(w*(0:m-1))*a);
t=w(i);