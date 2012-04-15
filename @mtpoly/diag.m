function q=diag(p)
%

% AM 09.01.09

n=max(p.n,p.m);
m=min(p.n,p.m);
s=p.s;
if m==1,
    d=max(s(:,1),s(:,2));
    s(:,1:2)=[d d];
    q=mtpoly(n,n,s);
else
    s=s(s(:,1)==s(:,2),:);
    s(:,2)=1;
    q=mtpoly(m,1,s);
end
