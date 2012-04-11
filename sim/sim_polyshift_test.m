function sim_polyshift_test

m=3+round(7*rand);
p=randn(1,m);
a=randn;
q=sim_polyshift(p,a);
t=randn;
fprintf(' sim_polyshift: m=%d, a=%f, t=%f, ',m,a,t)
fprintf('error=%f\n',abs(polyval(p,t+a)-polyval(q,t)))