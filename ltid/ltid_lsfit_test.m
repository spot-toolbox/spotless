function ltid_lsfit_test(m)

if nargin<1, m=1000; end

w=linspace(0.5,1,100)';
t=linspace(0,pi,5000)';
z=exp(1i*w);
z1=0.98*exp(1i*0.7);
z2=0.95*exp(1i*1.1);
z3=0.99*exp(1i*0.1);
u=z./(z-1)+(z/z1)./(z/z1-1)+(z/z1')./(z/z1'-1)+ ...
    (z/z3)./(z/z3-1)+(z/z3')./(z/z3'-1)+ ...
    (z/z2)./(z/z2-1)+(z/z2')./(z/z2'-1);

v=ltid_lsfit(u,w,m,t);

close(gcf);
subplot(2,1,1);plot(w,real(u),'.',t,real(v));grid
subplot(2,1,2);plot(w,imag(u),'.',t,imag(v));grid