function a=ltid_q2a(q)
% function a=ltid_q2a(q)
%
% given 1-by-(m+1) real q, find (m+1)-by-1 real a such that
% cos(t*(0:m))*a=q*exp((m:-1:0)'*1i*t) for all real t

if nargin<1, error('1 input required'); end
if ~isa(q,'double'), error('input 1 not a double'); end
a=real(q(:));
m=length(a)-1;
a=conv(a,a(m+1:-1:1));
a=[a(m+1);2*a(m+2:2*m+1)];