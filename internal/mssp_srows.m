function s=mssp_srows(s0)
% function s=mssp_srows(s0)
%
% order rows of an mssp structure matrix variable-wise

s=s0;
[ms,ns]=size(s);
if ns<7, return; end   % nothing to do if less than 2 variables per term
k=round((ns-3)/2);

[y,ii]=sort(s(:,3:2+k),2,'descend');     % order rows variable-wise
%[y,ii]=sort(-s(:,3:2+k),2);
sd=s(:,3+k:2+2*k);
ii=repmat((1:ms)',1,k)+ms*(ii-1);
s=[s(:,1:2) y sd(ii) s(:,ns)];
