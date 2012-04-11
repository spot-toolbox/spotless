function v=ltid_abw2v(a,b,w)
% function v=ltid_abw2v(a,b,w)
%
% for (m+1)-by-1 a, (m+1)-by-k b, and n-by-1 w,
% v=(cos(w*(0:m))*b)./repmat(cos(w*(0:m))*a,1,k)

m=size(a,1)-1;
k=size(b,2);
cs=cos(w(:)*(0:m));
v=(cs*b)./repmat(cs*a,1,k);
