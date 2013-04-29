function [n,nf,nl,nq,nr,ns] = spot_sdp_cone_dim(K)
    nf = 0;
    nl = 0;
    nq = 0;
    nr = 0;
    ns = 0;

    
    if ~isempty(K.f)
        nf = K.f;
    end

    if ~isempty(K.l)
        nl = K.l;
    end
    
    if ~isempty(K.q)
        nq = sum(K.q);
    end

    if ~isempty(K.r)
        nr = sum(K.r);
    end
    
    if ~isempty(K.s)
        ns = sum(K.s.^2);
    end

    n = nf + nl + nr + nq + ns;
    
end