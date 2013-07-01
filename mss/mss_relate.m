function ik=mss_relate(x,y)
% function ik=mss_relate(x,y)
%
% for two columns x,y the rows of the two-column matrix ik 
% stores all pairs (i,k) such that x(i)=y(k)
% works on data convertible to "double"


x=double(x);
y=double(y);
x=x(:);
y=y(:);
mx=length(x);
my=length(y);
if mx*my==0, ik=[]; return; end
m=mx+my;
g=sortrows([[x  -(1:mx)'];[y (1:my)']]);
g=[[g(1:m-1,1)==g(2:m,1);0] g(:,2)];

ik=spot_gset(g);