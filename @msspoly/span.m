function q = span(p)
%
% q = span(p)
%
%   p -- n-by-1 msspoly
%   q -- m-by-1 msspoly, m <= n
%
%  Polynomials {q(i)}_i have the same span as those {p(i)}_i
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
    
    qCoeff = Coeff;
    
    for i = 1:n
        % Term is not the zero poly.
        if ~all(Coeff(i,:)==0)
            if m == 0
                m = 1;
                qCoeff(1,:) = Coeff(i,:);
            else
                b = Coeff(i,:);
                A = qCoeff(1:m,:);
                ii = ~all([ A ; b] == 0);
                c = b(1,ii)/A(:,ii);
                if norm(b(1,ii)-c*A(:,ii))/norm(b) > 1e-6
                    m = m+1;
                    qCoeff(m,:) = Coeff(i,:);
                end
            end
        end
    end
    
    q = recomp(var,pow,qCoeff(1:m,:));
end