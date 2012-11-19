function [prog,S] = psd(prog,Q)
    if size(Q,1) ~= size(Q,2)
        error('Second argument must be square matrix.');
    end
    [prog,S] = new(prog,size(Q,2),'psd');
    prog = eq(prog,mss_s2v(S-Q));
end