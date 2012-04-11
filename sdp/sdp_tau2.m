function t=sdp_tau2(p)
% function t=sdp_tau2(p)
%
% For n-by-3 p (p(:,3)>0) yields min{t>0: sign(p*[t^2;t;1]) changes at t}

%if min(p(:,3))<=0, fprintf('\n min(p)=%e\n',min(p(:,3))); end
d=p(:,2).^2-4*(p(:,1).*p(:,3));
s=(d>max(0,p(:,2).*abs(p(:,2))));
t=min(2*p(s,3)./(sqrt(d(s))-p(s,2)));