function ltid_q2a_test(m)

if nargin<1, m=5; end
q=randn(m+1,1);
a=ltid_q2a(q);
t=linspace(0,pi,100);
cs=cos(t'*(0:m));
ex=exp(t'*(1i*(m:-1:0)));
er=norm(cs*a-abs(ex*q).^2);
fprintf(' error: %f\n',er)