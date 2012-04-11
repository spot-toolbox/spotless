function [fun,coeff] = nlid_generate_basis(var,l,mn,orth)
    if nargin < 4, orth = []; end
    if size(mn,1) ~= 1, error('Monomials must be 1xK'); end
    if strcmp(orth,'legendre')
        % Form Gram matrix:
        K = mn'*mn;
        Ktr = mss_s2v(K);
        [xx,pp,MM] = decomp(Ktr);
        n = length(xx);
        % Data is in [-1,1]^n
        c = ones(size(pp,1),1);
        for i = 1:n
            k = pp(:,i);
            c = c.*(1./(k+1)).*(2*(mod(k+1,2)==1));
        end
        K = mss_v2s(MM*c);
        G = chol(inv(K));
        mn = mn*G';
    end
    if size(mn,1) == 1, mn = repmat(mn,l,1); end
%     elseif size(mn,1) ~= l, error([var ' monomials size ' ...
%                             'incorrect.']); 
%     end
    coeff = reshape(msspoly(var,prod(size(mn))),size(mn,1),size(mn,2));
    fun = sum(coeff.*mn,2);
end