function v=nlid_tpls(x,y,d,e)
% function v=nlid_tpls(x,y,d,e)
%
% INPUTS:
%   x  -  m-by-n real
%   y  -  m-by-k real
%   d  -  1-by-n positive integer
%   e  -  even tp flag (default e=0)
%
% OUTPUT:
%   v  -  m-by-k real
%
% v is the best approximation of y, row-wise, 
% by a trigonometric polynomial with degrees bounded by d

if nargin<3, error('3 inputs required'); end
if nargin<4, e=0; end
x=real(double(x));
y=real(double(y));
d=max(1,round(real(double(d(1,:)))));
[m,n]=size(x);
if m~=size(y,1), error('inputs 1,2 incompatible'); end
if n~=length(d), error('inputs 1,3 incompatible'); end
U=sim_tprange(x,d,e);
v=U*(U'*y);