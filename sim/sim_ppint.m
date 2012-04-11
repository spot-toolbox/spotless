function q=sim_ppint(p)
% function q=sim_ppint(p)
%
% integrate MATLAB's piecewise polynomial vector

if nargin<1, error('one input required'); end
if ~sim_isppvec(p), error('input is not a pp vector'); end
q=p;
m=q.dim;
d=q.order-1;
t=q.breaks;
n=length(t)-1;
c=q.coefs.*repmat(1./(d+1:-1:1),n*m,1);
T=t(2:n+1)-t(1:n);
h=sum(c.*kron(repmat(T(:),1,d+1).^repmat(d+1:-1:1,n,1),ones(m,1)),2);
h=reshape(h,m,n);
h=cumsum([zeros(m,1) h(:,1:n-1)],2);
q.coefs=[c h(:)];
q.order=d+2;
