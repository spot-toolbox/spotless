function q=ltid_a2q(a)
% function q=ltid_a2q(a)
%
% given real column a, find real column q of the same length such that 
% cos((0:m)*t)*a=const|exp((m:-1:0)*j*t)*q|^2

a=a(:);
m=length(a)-1;
if a(m+1)==0, a(m+1)=1e-3; end
p=[a(m+1:-1:2);2*a(1);a(2:m+1)];
z=roots(p);
y=sort(abs(z));
z=z(abs(z)<=y(m));
z=z./max(1,1.0001*abs(z));
q=real(poly(z));
q=q*(sqrt(a(1))/norm(q));