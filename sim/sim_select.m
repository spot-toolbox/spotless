function [XX,YY]=sim_select(X,Y,K)
% function [XX,YY]=sim_select(X,Y,K)
%
% given n pairs (xi,yi) where xi are row m-vectors, yi are row kqvectors
% stored in X=[x1;x2; ... xn], y=[y1;y2; ... yn],
% and an m-vector of positive integers K=[k1 k2 ... km],
% select (up to) one sample per bin, in the partition of the
% hypercube min(X,[],1)<x<max(X,[],1) into k1*k2*...*km bins

[n,m]=size(X);
Xmin=min(X,[],1);
Xmax=max(X,[],1);
W=(X-repmat(Xmin,n,1))./repmat((Xmax-Xmin)./K,n,1);
Z=round(W+0.5)-0.5;

%close(gcf);plot(Z(:,1),Z(:,2),'.');
KK=prod(toeplitz(ones(m,1),[1 K(1:m-1)]));
V=Z*KK'-0.1*abs(W(:,1)-Z(:,1));
[V,I]=sort(V);
%close(gcf);plot(1:n,V,'.');
XX=X(I(1:n-1),:);
YY=Y(I(1:n-1),:);
S=(0.5+V(1:n-1)<V(2:n));
XX=XX(S,:);
YY=YY(S,:);