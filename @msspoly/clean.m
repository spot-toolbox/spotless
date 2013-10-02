function a=clean(a,tol)
   if nargin < 2, tol = 1e-6; end
   [x,p,M]=decomp(a);
   M(abs(M)<tol)=0;
   a=recomp(x,p,M,size(a));
end
