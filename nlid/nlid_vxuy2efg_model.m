function [prog,model] = nlid_vxuy2efg_model(prog,n,m,k,domain,...
                                                emonom,fmonom, ...
                                                gmonom,affine)
    if nargin < 5, error('Must provide domain.'); end
    if nargin < 6, fmonom = @(x,u) [x;u;1]; end
    if nargin < 7, emonom = @(x) x; end
    if nargin < 8, gmonom = @(x,u) [x;u;1]; end
    if nargin < 9
        affine = struct('X',struct('mu',zeros(n,1),'T',eye(n)),...
                        'V',struct('mu',zeros(n,1),'T',eye(n)),...
                        'Y',struct('mu',zeros(k,1),'T',eye(k)),...
                        'U',struct('mu',zeros(m,1),'T',eye(m)));
    end
    
    v = msspoly('v',n);
    x = msspoly('x',n);
    u = msspoly('u',m);
    y = msspoly('y',k);
    
    prog = mssprog;
    
    [f,fcoeff] = nlid_vxuy2efg_generate_basis('f',n,fmonom(x,u)','legendre');
    [e,ecoeff] = nlid_vxuy2efg_generate_basis('e',n,emonom(x)','legendre');
    [g,gcoeff] = nlid_vxuy2efg_generate_basis('g',k,gmonom(x,u)','legendre');
    
    prog.free = fcoeff;
    prog.free = ecoeff;
    prog.free = gcoeff;
    
    model = struct('domain',domain,'v',v,'x',x,'u',u,'y',y,...
                   'e',e,'f',f,'g',g,...
                   'affine',affine);
end
                                              