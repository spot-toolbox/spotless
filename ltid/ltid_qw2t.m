function t=ltid_qw2t(q,w,K)
% function t=ltid_qw2t(q,w,K)
%
% DT frequency samples associated with Schur polynomial q and original w

e=roots(q);
if any(abs(e)>=1), error('not a Schur polynomial'); end
t=repmat(angle(e),1,2*K+1)+sqrt(K)*(1-abs(e))*((-K:K)/K);
t=sort([linspace(0,pi,10*K)';t(:);w(:)]);
t=t((t>=0)&(t<=pi));
n=length(t);
t=t([-1;t(1:n-1)]<t);