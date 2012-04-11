function q=sim_polyshift(p,a)
% function q=sim_polyshift(p,a)
%
% p,q are polynomals such that q(t)=p(t+a) for all t

n=size(p,2);        % degree
i=repmat(n:-1:1,n,1);
j=repmat((n:-1:1)',1,n);
q=p*((gamma(j)./(gamma(i).*gamma(j-i+1))).*toeplitz([1;zeros(n-1,1)],a.^(0:1:n-1)));