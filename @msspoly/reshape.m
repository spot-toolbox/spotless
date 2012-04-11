function y=reshape(x,m,n)

mx=size(x,1);
[ms,ns]=size(x.s);
h=x.s(:,1)+mx*(x.s(:,2)-1);
ii=mod(h-1,m)+1;
jj=round((h-ii)/m)+1;
y=msspoly(m,n,[ii jj x.s(:,3:ns)]);
