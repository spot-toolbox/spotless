function [prog,model] = nlid_vxug2ef_model(prog,n,m,g,domain,...
                                             emonom,fmonom,affine)
    if nargin < 5, error('Must provide domain.'); end
    if nargin < 6, emonom = @(x,u) [x;u;1]; end
    if nargin < 7, fmonom = @(x) x; end
    if nargin < 8
        affine = struct('X',struct('mu',zeros(n,1),'T',eye(n)),...
                        'V',struct('mu',zeros(n,1),'T',eye(n)),...
                        'Y',struct('mu',zeros(k,1),'T',eye(k)),...
                        'U',struct('mu',zeros(m,1),'T',eye(m)));
    end
    

    
    v = msspoly('v',n);
    x = msspoly('x',n);
    u = msspoly('u',m);
    k = size(g(x,u),1);
    y = msspoly('y',k);
    
    prog = mssprog;
    
    [f,fcoeff] = nlid_vxuy2efg_generate_basis('f',n,fmonom(x,u)','legendre');
    [e,ecoeff] = nlid_vxuy2efg_generate_basis('e',n,emonom(x)','legendre');
    
    prog.free = fcoeff;
    prog.free = ecoeff;
    
    % Modify g to match affine transformation:
    g = inv(affine.Y.T)*(g(affine.X.T*x+affine.X.mu,affine.U.T*u+affine.U.mu) - affine.Y.mu);
    
    model = struct('domain',domain,'y',y,'v',v,'x',x,'u',u,...
                   'e',e,'f',f,'g',g,...
                   'ecoeff',ecoeff,'fcoeff',fcoeff,'gcoeff',[],...
                   'affine',affine);
end
                                              