n=2;
x = msspoly('x',2);
c = msspoly('c',n);
s = msspoly('s',n);
cs = [c s];
phi=monomials([c;s],0:4);
p=randn(1,length(phi))*phi*(1+x'*x);

pr = spotsosprog;
pr = pr.withIndeterminate(x);
pr = pr.withIndeterminate([c;s]);
[pr,r] = pr.newFree(1);

mn = mpmonomials([c;s],0:4,x,0:2);

[pr,si] = pr.newSOSTrigPoly(mn,cs,mn,n);


pr = withSOSTrigPoly(pr,p-r-c'*si,cs,mn);

sol = pr.minimize(-r,@spot_mosek);

N=1e4;
Th = (pi/2)*(2*(rand(2,N)-0.5));
X = randn(n,N);
[v,I]=min(msubs(p,[c;s;x],[cos(Th) ; sin(Th) ; X]));
full([double(sol.eval(r)) v double(sol.eval(r))<v])