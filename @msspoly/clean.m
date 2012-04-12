function a=clean(a,tol)

[x,p,M]=decomp(a);
M(abs(M)<tol)=0;
a=recomp(x,p,M);
