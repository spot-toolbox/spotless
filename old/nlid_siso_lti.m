function [G,c,e]=nlid_siso_lti(X,m,k)
% function [G,c,e]=nlid_siso_lti(X,m,k)
%
% Optimally fit a DT SISO model y=Gu+c, where G is the LTI system
%   y(t)+a1*y(t-1)+a2*y(t-2)+...+am*y(t-m)=b0*u(t)+b1*u(t-1)+...+bk*u(t-k),
% stability of which is certified by the passivity constraint
%     Re[1+a1/z+a2/z^2+...+am/z^m]>0 for all |z|=1,
% and c is a constant, to the i/o data given in 
% X={[U1 Y1],...,[Un Yn]}, Ui=[ui(1);...;ui(Ti)], Yi=[yi(1);...;yi(Ti)].
%
% INPUTS:
%   X  -  a cell array of Ti-by-2 real matrices
%   m  -  a non-negative integer
%   k  -  a non-negative integer 
% 
% OUTPUTS:
%   G  -  SISO LTI DT model G(z)=(b0+b1/z+...+bk/z^k)/(1+a1/z+...+am/z^m)
%   c  -  a real number
%   e  -  error upper bound

if nargin<3, error('3 inputs required'); end
if ~iscell(X), error('input 1 must be a cell array'); end
m=max(0,round(real(double(m(1)))));
k=max(1,round(real(double(k(1)))));
[vv,ww,emsg]=nlid_siso_c2m(X,m,k);
if ~isempty(emsg), error(emsg); end
[f,e]=nlid_siso(ww,vv);
z=[msspoly('y',m+1);msspoly('u',k+1)];
a=double(diff(f,z));
c=-double(subs(f,z,zeros(size(z))))/sum(a(1:m+1));
if m>k,
    G=-tf([a(m+2:m+2+k) zeros(1,m-k)],a(1:m+1),-1);
else
    G=-tf(a(m+2:m+2+k),[a(1:m+1) zeros(1,k-m)],-1);
end