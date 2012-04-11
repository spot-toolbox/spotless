function v=ltid_pqw2v(p,q,w)
% function v=ltid_pqw2v(p,q,w)
%
% for (m+1)-by-1 a, (m+1)-by-k b, and n-by-1 w,
% v=(zz*p')./repmat(zz*q',1,k), where zz=exp(w*j*(m:-1:0))

m=size(q,2)-1;
k=size(p,1);
zz=exp(w(:)*(1i*(m:-1:0)));
v=(zz*p')./repmat(zz*q',1,k);
