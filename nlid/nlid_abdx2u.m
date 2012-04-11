function u=nlid_abdx2u(a,b,d,x)
if (max(x(:))>1)||(min(x(:))<-1), 
    disp([min(x) max(x)])
    warning('argument not in [-1 1]'); 
end
x=max(min(x,1),-1);
cs=cos(acos(x)*d);
u=(cs*b)./(cs*a);