function y=sum(x,b)

if isempty(x), 
    y=x; 
    return; 
end
if nargin<2, 
    if size(x,1)==1,
        b=2;
    else
        b=1;
    end
end
s=x.s;
if b==1, 
    s(:,1)=1;
    y=mtpoly(1,x.n,s);
else
    s(:,2)=1; 
    y=mtpoly(x.m,1,s);
end
