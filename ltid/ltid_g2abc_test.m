function ltid_g2abc_test(n,m)

if nargin<1, n=100; end
if nargin<2, m=5; end
A=randn(m);
A=A/(1.1*max(abs(eig(A))));
B=randn(m,1);
C=randn(1,m);
g=impulse(ss(A,B,C,0,-1),n+1);
g=g(2:n+1);
g=g(:)+0.001*randn(n,1);

[a,b,c]=ltid_g2abc(g,m,'fire');
ge=impulse(ss(a,b,c,0,-1),n+1);
ge=ge(2:n+1);
close(gcf);plot(1:n,g,'.',1:n,ge); grid;pause

[a,b,c]=ltid_g2abc(g,m);
ge=impulse(ss(a,b,c,0,-1),n+1);
ge=ge(2:n+1);
close(gcf);plot(1:n,g,'.',1:n,ge); grid;pause

[a,b,c]=ltid_g2abc(g,m,'fir');
ge=impulse(ss(a,b,c,0,-1),n+1);
ge=ge(2:n+1);
close(gcf);plot(1:n,g,'.',1:n,ge); grid
