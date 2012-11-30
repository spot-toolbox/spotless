function q = span(p)
%
% q = span(p)
%
%   p -- n-by-1 msspoly
%   q -- m-by-1 msspoly, m <= n
%
%  Polynomials {q(k)}_k are a maximal independent subset of {p(i)}_i
%
    if size(p,2) ~=1
        error('Argument must be a column.');
    end
    
    n = size(p,1);
    q = msspoly(zeros(n,1));
    m = 0;
    
    [var,pow,Coeff] = decomp(p);
    
    % zero polynomials (or empty polynomials :\)
    if isempty(Coeff)
        q = [];
        return;
    end
    
    
    [Q,R] = qr(Coeff');
    
    Rsave = zeros(size(R));
    m = 1;
    for i = 1:n
        if R(m,i) ~= 0
            Rsave(:,m) = R(:,i);
            m = m+1; % next row!
        end
        if m > size(R,1), break; end
    end
    
    qCoeff = Rsave(:,1:m-1)'*Q';
    
    q = recomp(var,pow,qCoeff);
end