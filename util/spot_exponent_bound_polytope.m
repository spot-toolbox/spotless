function [ptope] = spot_exponent_bound_polytope(pow)
%
%   [ptope] = spot_exponent_bound_polytope(pow)
%
%   pow   -- matrix of nonnegative integers.
%   ptope -- matrix of nonnegative integers such that
%            ZZ^n intersect conv(pow(:,1),pow(:,2),...,pow(:,end))
%            belongs to ptope.
%
    hpow = unique(floor(pow/2),'rows');
    
    n = size(pow,2);
    d = max(sum(hpow,2));
    
    % Generate [0,d]^n bounding box.
    mpow = mss_asd(n+1,d);
    ptope = mpow(:,2:end);
    
    % Quickly remove some obviously unnecessary monomials (i.e.
    % total degree or scalar degree too high / low
    L = [ ones(n,1) eye(n) -eye(n) -ones(n,1)];
    opt = max(hpow*L,[],1);
    ptope = ptope(all(ptope*L <= repmat(opt,size(ptope,1),1),2),:);
end