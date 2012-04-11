function G=nlid_io_lti(X,m,k)
% function G=nlid_io_lti(X,m,k)
%
% Optimally fit a DT SISO LTI model (input u, output y) of the form
% a0*y(t)+a1*y(t-1)+a2*y(t-2)+...+am*y(t-m)=b1*u(t-1)+...+bk*u(t-k),
% stability of which is certified by the passivity constraint
%     Re[a0+a1/z+a2/z^2+...+am/z^m]>0 for all |z|=1,
% to the i/o data given in X={[U1 Y1],[U2 Y2],...,[Un Yn]}, where
% Ui=[ui(1);ui(2);...;ui(Ti)], Yi=[yi(1);yi(2);...;yi(Ti)].
%
% INPUTS:
%   X  -  a cell array of Ti-by-2 real matrices
%   m  -  a non-negative integer (default m=0)
%   k  -  a positive integer (default k=max(m,1))
% 
% OUTPUTS:
%   G  -  SISO LTI DT model G(z)=(b1/z+...+bk/z^k)/(a0+a1/z+...+am/z^m)

if nargin<1, error('one input required'); end
if ~iscell(X), error('1st input must be a cell array'); end
if nargin<2, m=0; end
if nargin<3, k=max(m,1); end
[vv,ww,emsg]=nlid_io_c2m(X,m,k);
if ~isempty(emsg), error(emsg); end
v=msspoly('v',m+1);
w=msspoly('w',k);
d=msspoly('d',m+1);
F=[v;w];
U=[v;w;d];
R=monomials(F,0:2);
if m>0,
    a=nlid_miso(vv,ww,F,R,U,d(2:m+1));
else
    a=nlid_miso(vv,ww,F,R,U);
end
G=-tf([0 a(m+2:m+1+k)'],[1 zeros(1,k)],-1)/tf(a(1:m+1)',[1 zeros(1,m)],-1);
%G=-tf(a(m+1:m+k)',a(1:m)',-1);
