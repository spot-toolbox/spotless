% Polynomial pieces to Piecewise Polynomial
%
% pp = sim_ps2pp(t,p,ts)
%
% Input:
%    t  -- scalar free msspoly.
%    p  -- n-by-N msspoly function of t.
%    ts -- 1-by-N+1 monotonically increasing break points.
%
% Output:
%    pp -- piecewise polynomial s.t. 
%          pp(t) = { p(:,i)(t-ts(i))  t \in [ts(i),ts(i+1))
function pp = sim_ps2pp(t,p,ts)
    n = size(p,1);
    N = size(p,2);
    if size(ts,2)-1 ~= N, error('length of p, ts disagree'); end
    
    [x,e,C]=decomp(p);

    if size(x) == [0 1]
        x = t; e = 0;
        if size(C) == [0 1]
            C = zeros(N,1);
        elseif size(C) == [1 N]
            C = C';
        end
    end
        
    if ~(all(size(x) == [1 1]) &...
         double(x-t) == 0)
        error('p must be a function only of t');
    end
    o = length(e); % + 1?
    [es,I] = sort(e,'descend');

    pp.breaks = ts;
    pp.coefs = full(C(:,I));
    pp.order = o;
    pp.pieces = N;
    pp.dim = [n];
    pp.form = 'pp';
end