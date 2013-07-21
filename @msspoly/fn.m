function f = fn(p,x)
% fn -- Convert an msspoly to a function of doubles in a
%        specified order.
%
% f = fn(p,x)
%  
% Arguments:
%     p -- A polynomial in x.
%     x -- A free n-by-1 msspoly of unique degree 1 monomials.
%

    [xf,pf,Mf]=decomp(msspoly(p));
    
    if isempty(Mf)
        f = @(x) zeros(size(p));
        return;
    end

    if isempty(xf)
        f = @(x) reshape(Mf,size(p));
        return;
    end
    
    
    [~,xfid] = isfree(xf);
    [fr,xid]  = isfree(x);
    
    if ~fr, error('Second argument must be free'); end
    
    mtch = mss_match(xid,xfid);
    
    if any(mtch == 0)
        error(['First argument must be a function of elts. of second ' ...
               'argument']);
    end
    
    
    function r = preal(xreal)
        if isa(xreal,'msspoly')
            xx = indexinto(xreal,mtch);
        else
            xx = xreal(mtch);
        end
        r = reshape(Mf*prod(repmat(xx.',size(pf,1),1).^pf,2),size(p));
    end
    f = @preal;
end