function U=sim_tprange(x,d,e)
% function U=sim_tprange(x,d,e)
%
% INPUTS:
%   x  -  m-by-n real
%   d  -  1-by-n positive integer
%   e  -  even tp flag (default e=0)
%
% OUTPUT:
%   U  -  m-by-k real
%
% columns of U form an orthonormal basis in the range of
% trigonometric polynomials of degrees bounded by d applied to x row-wise

if nargin<2, error('2 inputs required'); end
if nargin<3, e=0; end
x=real(double(x));
d=max(1,round(real(double(d(1,:)))));
[m,n]=size(x);
if n~=length(d), error('inputs 1,3 incompatible'); end
x0=min(x,[],1);
x1=max(x,[],1);
dx=x1-x0+(x1==x0);
u=acos((x-repmat(x0,m,1))./repmat(dx,m,1));
dd=mint_ch(mint_down(d));
dd=dd(2:size(dd,1),:);
M=exp(u*(1i*dd'));
if e==1,
    U=orth([real(M) imag(M)]);
else
    U=orth(real(M));
end