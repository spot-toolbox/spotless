function G=ltid_reshape(G0,m,n)
% function G=ltid_reshape(G0,m,n)
%
% reshape for LTI systems

if nargin<3, error('3 inputs required'); end
if ~isa(G0,'lti'), error('input 1 not LTI'); end
m=max(1,round(real(m(1))));
n=max(1,round(real(n(1))));
[mg,ng]=size(G0);
if m*n~=mg*ng, error('incompatible dimensions'); end
[p,q]=tfdata(G0);
G1=tf(p(:),q{1,1},-1);
G=ss([],[],[],zeros(m,n),-1);
for i=1:n, G(:,i)=G1(1+m*(i-1):m*i); end

