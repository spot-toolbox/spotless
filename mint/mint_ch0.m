function y=mint_ch0(x,F,R)
% function y=mint_ch0(x,F,R)
%
% for a double integer matrices x(m-by-n)), F((n+1)-by-p), R(1-by-1)
% such that [1 x] is left invertible and R>0, 
% y lists all integer vectors v in the convex hull of the rows of x 
% such that mod([v 1]*F,R)=0

if nargin<1, error('one input required'); end
if ~mint_isint(x), error('1st input is not a double integer'); end
[m,n]=size(x);    
if rank([ones(m,1) x])<=n, error('singular input'); end
if nargin<3, F=zeros(n+1,1); R=1; end
if ~mint_isint(F), error('input 2 is not a double integer'); end
if ~mint_isint(R), error('input 3 is not a double integer'); end
if size(F,1)~=n+1, error('inputs 1,2 incompatible'); end
if ~isscalar(R), error('input 3 not a scalar'); end
if R<1, error('input 3 not positive'); end

if n==1, 
    y=(min(x):max(x))'; 
else
    f=mint_eqr(F,R);
    y0=mint_ch0(x(:,2:n),f,R); % y(t,:)=[r(t) y0(s(t),:)]
    k0=size(y0,1);
    [a,b]=mint_v2i(x);
    c0=zeros(k0,1);                       % c0(s(i))>=0 is necessary
    ee=(a(1,:)==0);
    if any(ee),   
        a0=a(:,ee);
        b0=b(ee);
        c0=min(repmat(b0,k0,1)-y0*a0(2:n,:),[],2);  
    end
    c1=repmat(Inf,k0,1);
    ee=(a(1,:)>0);
    if any(ee),                           % r(t)<=c1(s(t))
        a1=a(:,ee);
        b1=b(ee);
        c1=min(floor((repmat(b1,k0,1)-y0*a1(2:n,:))./repmat(a1(1,:),k0,1)),[],2);
    else
        save bugv2i x
        error('>')
    end
    c2=repmat(-Inf,k0,1);                 % r(t)>=c2(s(t))
    ee=(a(1,:)<0);
    if any(ee),
        a2=a(:,ee);
        b2=b(ee);
        c2=max(ceil((repmat(b2,k0,1)-y0*a2(2:n,:))./repmat(a2(1,:),k0,1)),[],2); 
    else
        save bugv2i x
        error('<')
    end
    e=(c0>=0)&(c1>=c2);
    k=sum(c1(e)-c2(e)+1);
    y=zeros(k,n);
    t=0;                        % number of finished rows in y 
    for s=1:k0,
        d=c1(s)-c2(s)+1;
        if (c0(s)>=0)&&(d>0),
            y(t+1:t+d,:)=[(c2(s):c1(s))' repmat(y0(s,:),d,1)];
            t=t+d;
        end
    end
end
ee=all(mod([y ones(size(y,1),1)]*F,R)==0,2);
y=y(ee,:);

