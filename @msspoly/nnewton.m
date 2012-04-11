function x=nnewton(p,v,k,x0,q)
% function x=nnewton(p,v,k,x0,q)
%
% for mss polynomial p in free msspoly variable v,
% performs k (default k=1) iterations 
%   x -> y=x-inv(q(x))*p(x) -> xnew=y./max(1,abs(y))
% starting at x=x0 

x=x0;
for i=1:k,
    x=x-inv(double(subs(q,v,x)))*double(subs(p,v,x));
    x=x./max(1,abs(x));
end