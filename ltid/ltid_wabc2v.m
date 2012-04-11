function v=ltid_wabc2v(w,a,b,c)
% function v=ltid_wabc2v(w,a,b,c)
%
% for n-by-1 w, (m+1)-by-1 a, (m+1)-by-k b, and (optional) m-by-k c
% v=(cos(w*(0:m))*b)./repmat(cos(w*(0:m))*a,1,k)

if nargin<4, c=[]; end
if isempty(c),
    m=size(a,1)-1;
    k=size(b,2);
    cs=cos(w(:)*(0:m));
    v=(cs*b)./repmat(cs*a,1,k);
elseif isempty(b),
    [m,k]=size(c);
    v=(1i*(sin(w(:)*(1:m))*c))./repmat(cos(w(:)*(0:m))*a,1,k);
else
    [m,k]=size(c);
    cs=cos(w(:)*(0:m));
    v=(cs*b+1i*(sin(w(:)*(1:m))*c))./repmat(cs*a,1,k);
end
