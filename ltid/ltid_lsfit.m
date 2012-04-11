function v=ltid_lsfit(u,w,m,t)
% function v=ltid_lsfit(u,w,m,t)
%
% least squares extrapolation
%
% INPUTS:
%   u:   n-by-k complex (frequency response)
%   w:   n-by-1 real (DT frequencies)
%   m:   positive integer (length of the frequency response)
%   t:   r-by-1 real (new DT frequencies, default t=w)
%
% OUTPUTS:
%   v:   r-by-k complex (extrapolated values)

if nargin<3, error('3 inputs required'); end
[n,k]=size(u);
if m<2*n+1, error('input 3 value is too small'); end
w=w(:);
if nargin<4, t=w; end
t=t(:);
r=length(t);

M=exp((-1i*w)*(0:m));
x=[real(M);imag(M);diag(0.1*(0:m))]\[real(u);imag(u);zeros(m+1,k)];
v=exp((-1i*t)*(0:m))*x;
