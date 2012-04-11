function q=sim_ppdiff(p)
% function q=sim_ppdiff(p)
%
% differentiate MATLAB's piecewise polynomial

if nargin<1, error('one input required'); end
if ~sim_isppvec(p), error('input is not a pp vector'); end
q=p;
if q.order<2, q.coefs=zeros(size(q.coefs)); return; end
q.order=q.order-1;
q.coefs=q.coefs(:,1:q.order).*repmat(q.order:-1:1,size(q.coefs,1),1);
