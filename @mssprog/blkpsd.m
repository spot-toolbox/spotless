function [prog,S] = blkpsd(prog,Q)
    n = size(Q,1);
    m=round((sqrt(1+8*n)-1)/2);
    if m*(m+1)~=2*n, error('Argument Two: wrong dimension for the psd type'); end
    
    [prog,S] = new(prog,[m size(Q,2)],'psd');
    prog = eq(prog,S-Q);
end