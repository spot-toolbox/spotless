function q=repmat(p,m,n)
%

% AM 09.01.09

[mp,np]=size(p);
mq=m*mp;
nq=n*np;
if isempty(p.s),
    q=msspoly(mq,nq,p.s);
    return
end
[ms,ns]=size(p.s);
s=repmat(p.s,m*n,1);
s(:,1)=s(:,1)+kron((0:m-1)',repmat(mp,n*ms,1));
s(:,2)=s(:,2)+repmat(kron((0:n-1)',repmat(np,ms,1)),m,1);
q=msspoly(mq,nq,s);



