function p=mtimes(p1,p2)
%

% AM 17.02.09

p1=msspoly(p1);
p2=msspoly(p2);
if (p1.n==1)&&(p1.m==1),
    p1=diag(repmat(p1,p2.m,1));
end
if (p2.n==1)&&(p2.m==1),
    p2=diag(repmat(p2,p1.n,1));
end
if p1.n~=p2.m, error('incompatible dimensions'); end
[m1,n1]=size(p1.s);
k1=round((n1-3)/2);
[m2,n2]=size(p2.s);
k2=round((n2-3)/2);
% new version begins
ik=mss_relate(p1.s(:,2),p2.s(:,1));
if isempty(ik),
    s=zeros(0,3);
else
    x=ik(:,1);
    y=ik(:,2);
    s1=p1.s;
    s2=p2.s;
    s=[s1(x,1) s2(y,2) s1(x,3:2+k1) s2(y,3:2+k2) ...
        s1(x,3+k1:2+2*k1) s2(y,3+k2:2+2*k2) s1(x,n1).*s2(y,n2)];
end
% old version begins
%s=[repmat(p1.s,m2,1) kron(p2.s,ones(m1,1))];
%s=s(s(:,2)==s(:,n1+1),:);
%s=[s(:,1) s(:,n1+2) s(:,3:2+k1) s(:,n1+3:n1+2+k2) ...
%    s(:,3+k1:2+2*k1) s(:,n1+3+k2:n1+2+2*k2) s(:,n1).*s(:,n1+n2)];
p=msspoly(p1.m,p2.n,s);